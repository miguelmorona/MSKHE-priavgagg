// Copyright 2026 Miguel Morona Mínguez
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.package main

package main

import (
	"flag"
	"log"
	"math"
	"math/big"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

// ----------------------------------------------------------------------------------------------//

// To check no errors happened during execution of a function
func check(err error) {
	if err != nil {
		panic(err)
	}
}

// To measure time execution of operations done by one party
func runTimed(f func()) time.Duration {
	start := time.Now()
	f()
	return time.Since(start)
}

// To measure time execution of the same operations which are simultaneously done in "parallel" by N parties
func runTimedParty(f func(), N int) time.Duration {
	start := time.Now()
	f()
	return time.Duration(time.Since(start).Nanoseconds() / int64(N))
}

// Definition of command line variables
var flagCommandLine = flag.Bool("commandline", false, "run the example with the command-line parameters.")
var flagNumParties = flag.Int("numparties", 20, "number of input parties.")
var flagNumCtxPerParty = flag.Int("numctxperparty", 10, "size of the model updates to aggregate.")
var flagQLevel = flag.Int("qlevel", 2, "level of Q.")
var flagPLevel = flag.Int("plevel", 0, "level of P.")
var flaglogN = flag.Int("logN", 10, "logarithm of lattice dimension.")
var flagBitLevelSize = flag.Int("bitlevelsize", 16, "size in bits of each level.")
var flagPreComputeA = flag.Bool("precomputea", false, "'a' polynomials components are precomputed before the encryption phase.")

type parameters struct {

	// Protocol parameters
	NumParties  int
	n           int // n = Number of ciphertexts per party to aggregate in each round (n*N corresponds to the "Size" of the model to aggregate)
	PreComputeA bool

	// Cryptographic parameters
	logN        int    // quotient polynomial Ring degree
	logQ        [2]int // logQ[0] = #Primes, logQ[1] = Primes bit-size
	plevel      int    // Maximum level for the modulus "p" (level 0 is the lowest available level)
	bench_Delta int    // Delta = 2^bench_Delta. It must be less than q1 (second limb) such that m~q0 and Delta*m ~ q1*q0 = plevel
	bench_eps   int    // Epsilon = 2^bench_eps. It must be less than q1 (second limb) such that m~q0 and Delta*m ~ q1*q0 = plevel
}

var benchParameters = []parameters{
	//----------------------------------------------------------------------------------------------------------------------//
	{NumParties: 32, n: 8000, PreComputeA: false, logN: 11, logQ: [2]int{2, 21}, plevel: 0, bench_Delta: 40, bench_eps: 30},
	//----------------------------------------------------------------------------------------------------------------------//
}

//-------------------------------------------------------------------------------------------------------------------------------//

//Set 1: --
//Set 2: {NumParties: 32, n: 8000, PreComputeA: false, logN: 11, logQ: [2]int{2, 21}, plevel: 0, bench_Delta: 40, bench_eps: 30},
//Set 3: --
//Set 4: {NumParties: 32, n: 4000, PreComputeA: false, logN: 12, logQ: [2]int{2, 46}, plevel: 0, bench_Delta: 90, bench_eps: 80},

//-------------------------------------------------------------------------------------------------------------------------------//

// Definition of variables used to measure runtime for specific steps of the aggregation protocol
var elapsedSetupCloud time.Duration // elapsedSetupCloud = 0 (no setup required for the Aggregator)
var elapsedSetupParty time.Duration
var elapsedInputParty time.Duration
var elapsedSetupParty_ri time.Duration
var elapsedSetupGaussianParty time.Duration
var elapsedPreprocessingEncryptParty time.Duration
var elapsedPreprocessingUniformSamplersParty time.Duration
var elapsedEncryptParty time.Duration
var elapsedEncryptPartyCPU time.Duration
var elapsedEncryptPartyWall time.Duration
var elapsedEncryptCloud time.Duration
var elapsedEvalCloudCPU time.Duration
var elapsedEvalCloud time.Duration
var elapsedEvalCloudWall time.Duration
var elapsedEvalParty time.Duration
var elapsedDecCloud time.Duration
var elapsedDecParty time.Duration
var elapsedDecCloudCPU time.Duration
var elapsedDecCloudWall time.Duration

var SetupMean time.Duration
var GenInputMean time.Duration
var EncPhaseMean time.Duration
var AggMean time.Duration
var DecMean time.Duration
var TotalMean time.Duration

// party: defines the inputs and private information of party P_i
type party struct {
	sk ring.Poly // sk_i (individual secret key of party P_i)
	r  ring.Poly // r_i (individual secret share of party P_i)

	input                    []ring.Poly
	PaddedNumModelParameters int // Size of the vector of model updates with zero padding (till the next multiple of params.N() => PaddedNumModelParameters = k*params.N(), where k is the smallest integer satisfying k*params.N() >= real number of NumModelParameters)
}

// AggRings: structure used to store variables related to the used quotient rings in the aggegation protocol
type AggRings struct {
	ringQ *ring.Ring // Z_{Q}/(X^N+1)
	Delta *big.Int   // Δ
}

// newAggRings: defines and initializes an AggRings variable
func newAggRings(params parameters) *AggRings {

	var err error
	N := 1 << params.logN // N is set as 2^logN
	rings := new(AggRings)

	// Create a generator for primes (of the specified bit size) compatible with negacyclic NTT for the given N
	g := ring.NewNTTFriendlyPrimesGenerator(uint64(params.logQ[1]), uint64(2*N))

	// Generate k = params.logQ[0] NTT-friendly primes, each approximately close to 2^logQ[1]
	primes, err := g.NextAlternatingPrimes(params.logQ[0])
	check(err)

	// Create the polynomial ring Z[x]_Q / (x^N + 1) with a modulus over an RNS with as many primes as in "primes"
	rings.ringQ, err = ring.NewRing(N, primes)
	check(err)

	rings.Delta = new(big.Int).Lsh(big.NewInt(1), uint(params.bench_Delta)) // Param set 2. Delta = 2^40. Param set 3. Delta  = 2^90

	return rings
}

// lowNormSampler: this structure is used to store the gaussian values in *Big.Int. We indicate the ring base in which we work (mod Q =  limb0*limb1*...)
type lowNormSampler struct {
	baseRing *ring.Ring
	coeffs   []*big.Int
}

// newLowNormSampler: initializes a lowNormSampler structure based on the provided base ring. Given a "baseRing" of type *ring.Ring, it allocates memory for the lowNormSampler structure.
func newLowNormSampler(baseRing *ring.Ring) (lns *lowNormSampler) {
	lns = new(lowNormSampler)
	lns.baseRing = baseRing
	lns.coeffs = make([]*big.Int, baseRing.N())
	return
}

// newPolyLowNorm genera un polinomio con coeficientes aleatorios uniformes en el rango cerrado [-delta, delta].
func (lns *lowNormSampler) newPolyLowNorm(delta *big.Int) (pol ring.Poly) {

    pol = lns.baseRing.NewPoly()
    prng, _ := sampling.NewPRNG()

    rangeSize := new(big.Int).Mul(delta, big.NewInt(2))
    rangeSize.Add(rangeSize, big.NewInt(1))

    for i := range lns.coeffs {
        tmp := bignum.RandInt(prng, rangeSize)        
        lns.coeffs[i] = new(big.Int).Sub(tmp, delta)
    }

    lns.baseRing.AtLevel(pol.Level()).SetCoefficientsBigint(lns.coeffs, pol)
    return
}

func PolyToFloat64Slice(pol ring.Poly, r *ring.Ring) []float64 {
	level := pol.Level()
	N := r.N()
	bigCoeffs := make([]*big.Int, N)
	coeffs := make([]float64, N)
	r.AtLevel(level).PolyToBigint(pol, 1, bigCoeffs)
	for i := 0; i < N; i++ {
		f, _ := new(big.Float).SetInt(bigCoeffs[i]).Float64()
		coeffs[i] = f
	}
	return coeffs
}

// main function. Example execution on terminal:

// go run ./CKKS-MSK.go 

func main() {

	file, err := os.Create("ckks-results.txt")
	check(err)
	defer file.Close() 
	l := log.New(file, "", 0)

	//l := log.New(os.Stderr, "", 0)

	flag.Parse()
	paramsSets := benchParameters
	if *flagCommandLine {
		paramsSets = benchParameters[0:1]
		paramsSets[0].logN = *flaglogN
		paramsSets[0].logQ = [2]int{*flagQLevel + 1, *flagBitLevelSize}
		paramsSets[0].plevel = *flagPLevel
		paramsSets[0].n = *flagNumCtxPerParty
		paramsSets[0].NumParties = *flagNumParties
		paramsSets[0].PreComputeA = *flagPreComputeA
	}


	// ----------------------------------------------------------------------------------------------//
	//                        Running several aggregation examples
	// ----------------------------------------------------------------------------------------------//

	for _, param := range paramsSets {

		// ----------------------------------------------------------------------------------------------//
		//                   Cryptographic parameters and initialization of samplers
		// ----------------------------------------------------------------------------------------------//

		// Initialize cryptographic parameters and structures for aggregation
		priaggrings := newAggRings(param) // Aggregation ring setup
		ringQ := priaggrings.ringQ

		// Extract protocol and cryptographic parameters from the current configuration
		n := param.n // Number of ciphertexts per party to aggregate
		NumParties := param.NumParties   // Number of parties in the protocol
		PreComputeA := param.PreComputeA // Flag for pre-computation optimization
		_ = PreComputeA
		eps := param.bench_eps // Epsilon for precision in result comparison

		l.Printf("\n\t\t\t||=======================================================================================================================||")
		l.Printf("\t\t\t||                    PUSHING FORWARD MULTI-SECRET-KEY HOMOMORPHIC ENCRYPTION FOR PRIVATE AGGREGATION                    ||")
		l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
		l.Printf("\t\t\t||                                                    ▄▖▖▖▖▖▄▖                                                           ||")
		l.Printf("\t\t\t||                                                    ▌ ▙▘▙▘▚                                                            ||")
		l.Printf("\t\t\t||                                                    ▙▖▌▌▌▌▄▌                                                           ||")
		l.Printf("\t\t\t||                                                                                                                       ||")
		l.Printf("\t\t\t||                                               ▄▖▄▖▄▖▄▖▖  ▖▄▖▄▖▄▖▄                                                     ||")
		l.Printf("\t\t\t||                                               ▌▌▙▌▐ ▐ ▛▖▞▌▐ ▗▘▙▖▌▌                                                    ||")
		l.Printf("\t\t\t||                                               ▙▌▌ ▐ ▟▖▌▝ ▌▟▖▙▖▙▖▙▘                                                    ||") //miniwi
		l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
		l.Printf("\t\t\t||                                                                                                                       ||")
		l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
		l.Printf("\t\t\t||   1.a. SETUP                                                                                                          ||")
		l.Printf("\t\t\t||           - Number of parties involved in the process:  %2d                                                            ||", NumParties)
		l.Printf("\t\t\t||           - Dimension of the model:                     %8d                                                      ||", n*(1<<param.logN))
		l.Printf("\t\t\t||                                                                                                                       ||")
		l.Printf("\t\t\t||   1.b. CRYPTOGRAPHIC PARAMETERS                                                                                       ||")
		l.Printf("\t\t\t||           - Level of security:                          128 bits                                                      ||")
		l.Printf("\t\t\t||           - Ring degree:                                2^%2d                                                          ||", param.logN)
		l.Printf("\t\t\t||           - Ciphertext modulus (Q):                     %2d bits                                                       ||", uint(param.logQ[0]*param.logQ[1]))
		l.Printf("\t\t\t||           - Number of values per ciphertext (MaxSlots): %4d                                                          ||", (1 << (param.logN - 1)))
		l.Printf("\t\t\t||                                                                                                                       ||")

		SetupMean = time.Duration(0)
		GenInputMean = time.Duration(0)
		EncPhaseMean = time.Duration(0)
		AggMean = time.Duration(0)
		DecMean = time.Duration(0)

		nexperiments := 50

		for z := 0; z < nexperiments; z++ {

			l.Printf("Number of experiment: n=%2d", z)

			// Initialize a pseudorandom number generator for sampling
			prng, err := sampling.NewPRNG()
			check(err)

			// Initialize ternary sampler for ringQ, with P:2.0/3.0 which gives a ternary uniform distribution {-1, 0, 1} -> 1/3, 1/3, 1/3
			ternarySamplerMontgomeryQ, err := ring.NewSampler(prng, ringQ, ring.Ternary{P: 2.0 / 3.0}, true)
			check(err)

			// Initialize Gaussian sampler with a standard deviation of 3.2, upper bound set to approximately 6 * sigma
			gaussianSamplerQ, err := ring.NewSampler(prng, ringQ, ring.DiscreteGaussian{Sigma: 3.2, Bound: 19.2}, false) //B_sigma = 6*Sigma
			check(err)
			_ = gaussianSamplerQ // Placeholder to prevent unused variable error

			// Initialize uniform sampler for ringQ for uniformly random coefficients
			uniformSamplerQ, err := ring.NewSampler(prng, ringQ, ring.Uniform{}, false)
			check(err)

			// Initialize low-norm sampler for generating low-norm polynomials in ringQ
			lowNormUniformQ := newLowNormSampler(ringQ)

			l.Printf("\t\t\t||   1.c. GENERATION OF INPUT PARTIES AND KEYS:                                                                          ||")

			elapsedSetupParty = time.Duration(0)
			elapsedInputParty = time.Duration(0)
			elapsedSetupParty_ri := time.Duration(0)

			P := genParties(priaggrings, ternarySamplerMontgomeryQ, NumParties)
			l.Printf("\t\t\t||           - Generation of secret keys sk_i                   %12s                                             ||", elapsedSetupParty)

			elapsedSetupParty_ri += runTimed(func() {
				genSetupShare(priaggrings, uniformSamplerQ, P)
			})
			l.Printf("\t\t\t||           - Generation of secret randomness r_i              %12s                                             ||", time.Duration(2)*elapsedSetupParty_ri)
			l.Printf("\t\t\t||                                                                                                                       ||")
			l.Printf("\t\t\t||           ✅ Done. Time per party: %12s                                                                       ||", time.Duration(2)*elapsedSetupParty_ri+elapsedSetupParty)
			l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
			l.Printf("\t\t\t||   2. GENERATION OF INPUTS FOR VERIFICATION                                                                            ||")
			aggexparray := genInputs(priaggrings, lowNormUniformQ, n, P, param) // The expected result is the aggregation of all the party models (aggexp = m_0 + m_1 + ... + m_{L-1})
			l.Printf("\t\t\t||                                                                                                                       ||")
			l.Printf("\t\t\t||           ✅ Done. Time per party: %12s                                                                       ||", elapsedInputParty)

			l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
			l.Printf("\t\t\t||   3. ENCRYPTION PHASE                                                                                                 ||")

			// encPhase:
			// 	Inputs => "a", "m", "sk", "r"
			// 	Outputs => enc(a, m, "mask"), partialdec(a, sk)

			// Encrypt: a*(si + ri) + e + Δ*m[i]
			encInput, partialDec := encPhase(priaggrings, gaussianSamplerQ, uniformSamplerQ, n, P)
			l.Printf("\t\t\t||                                                                                                                       ||")
			l.Printf("\t\t\t||           ✅ Done. Time per party: %12s                                                                       ||", elapsedEncryptCloud+elapsedEncryptParty)

			//PENDING UPDATE receiving only inputs from P[i], uniformSamplerQ from the same crs and a fresh gaussianSamplerQ
			//PENDING UPDATE run loop "for i := 0; i < n; i++ { encInputi, particDeci := encPhase(priaggrings, gaussianSamplerQ, uniformSamplerQ, n, P[i], param) }"

			l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
			l.Printf("\t\t\t||   4. AGGREGATION PHASE                                                                                                ||")

			// Aggregation phase:
			// Inputs => c0, c1, ..., c_{L - 1} (ciphertexts from all parties)
			// Outputs => cagg = c0 + c1 + ... + c_{L - 1} (aggregated ciphertext result)

			encShareAgg := evalPhase(priaggrings, n, encInput)
			l.Printf("\t\t\t||                                                                                                                       ||")
			l.Printf("\t\t\t||           ✅ Done. Total time: %12s                                                                           ||", elapsedEvalCloud+elapsedEvalParty)
			l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
			l.Printf("\t\t\t||   5. DECRYPTION PHASE                                                                                                 ||")

			//decPhase:
			// Inputs => partialdec(a, sk_i) (partial decryption of each party's input), aggoutputenc (aggregated encrypted ciphertext)
			// Outputs => recAggShare = partialdec(a,sk_i) + aggoutputenc (reconstructed aggregated result after decryption)º

			recAggShare := decPhase(priaggrings, partialDec, encShareAgg, n, param, P)

			var flatAggExp []float64
			for _, poly := range aggexparray {
				flatAggExp = append(flatAggExp, PolyToSignedFloat64Slice(poly, ringQ)...)
			}

			var flatRecAggShare []float64

			elapsedDecCloud += runTimedParty(func() {			
			for _, poly := range recAggShare {
				flatRecAggShare = append(flatRecAggShare, PolyToSignedFloat64Slice(poly, ringQ)...)
			}

			},len(P))

			l.Printf("\t\t\t||                                                                                                                       ||")
			l.Printf("\t\t\t||           ✅ Done. Time per party: %12s                                                                       ||", elapsedDecCloud+elapsedDecParty)
			l.Printf("\t\t\t||-----------------------------------------------------------------------------------------------------------------------||")
			l.Printf("\t\t\t||   6. PRECISION ANALYSIS                                                                                               ||")
			l.Printf("\t\t\t||           - Total running time (per client + agg): %12s                                                       ||", elapsedSetupParty_ri+elapsedSetupParty+elapsedEncryptCloud+elapsedEncryptParty+elapsedEvalCloud+elapsedEvalParty+elapsedDecCloud+elapsedDecParty)
			l.Printf("\t\t\t||           - Results decoded. Showing some rows:                                                                       ||")
			l.Printf("\t\t\t||                                                                                                                       ||")

			partiesBig := big.NewInt(int64(NumParties))
			modulus := new(big.Int).Mul(priaggrings.Delta, partiesBig)
			
			flatAggExp = DivideByBigInt(flatAggExp, modulus)
			flatRecAggShare = DivideByBigInt(flatRecAggShare, modulus)				

			maxAbsErr := 0.0
			nerror := 0
			relErrorThreshold := math.Exp2(float64(-eps))


			for i := range flatAggExp {
				diff := math.Abs(float64(flatAggExp[i] - flatRecAggShare[i]))

				// Track maximum error
				if diff > maxAbsErr {
					maxAbsErr = diff
				}

				if diff > relErrorThreshold {
					nerror++
				}			

				if i < 65 || i >= (len(flatAggExp))-65 {
					l.Printf("\t\t\t||   Row [%8d] of Aggregated Model w.  Exp Agg Result: %14.8f, Decrypted: %14.8f, AbsErr= %.2e  ||", i+1, flatAggExp[i], flatRecAggShare[i], diff)
				}
				

			}

			// If not all results match, log fail.
			if maxAbsErr > relErrorThreshold {
				l.Printf("\t\t\t||                                                                                      ||===============================||")
				l.Printf("\t\t\t||                                                                                      ||  Max Absolute Error: %.2e ||", maxAbsErr)
				l.Printf("\t\t\t||                                                                                      || Precision Threshold: %.2e ||", relErrorThreshold)
				l.Printf("\t\t\t||                                                                                      || Number of errors:    %7d  ||", nerror)
				l.Printf("\t\t\t||                                                                                      ||===============================||")
				l.Printf("\t\t\t||=======================================================================================================================||")
				l.Printf("\t\t\t||                                            ❌ ❌  INCORRECT :((                                                       ||")
				l.Printf("\t\t\t||=======================================================================================================================||")				
			} else {
				// If all results match, log success.
				l.Printf("\t\t\t||                                                                                      ||===============================||")
				l.Printf("\t\t\t||                                                                                      ||  Max Absolute Error: %.2e ||", maxAbsErr)
				l.Printf("\t\t\t||                                                                                      || Precision Threshold: %.2e ||", relErrorThreshold)
				l.Printf("\t\t\t||                                                                                      || Number of errors:    %7d  ||", nerror)
				l.Printf("\t\t\t||                                                                                      ||===============================||")
				l.Printf("\t\t\t||=======================================================================================================================||")
				l.Printf("\t\t\t||                                                  ✅✅ CORRECT!!                                                       ||")
				l.Printf("\t\t\t||=======================================================================================================================||")
			}

			

			// ----------------------------------------------------------------------------------------------//

			// ----------------------------------------------------------------------------------------------//

			SetupMean += (time.Duration(2)*elapsedSetupParty_ri + elapsedSetupParty) / time.Duration(nexperiments)
			GenInputMean += (elapsedInputParty) / time.Duration(nexperiments)
			EncPhaseMean += (elapsedEncryptCloud + elapsedEncryptParty) / time.Duration(nexperiments)
			AggMean += (elapsedEvalCloud + elapsedEvalParty) / time.Duration(nexperiments)
			DecMean += (elapsedDecCloud + elapsedDecParty) / time.Duration(nexperiments)
			TotalMean += (time.Duration(2)*elapsedSetupParty_ri + elapsedSetupParty + elapsedEncryptCloud + elapsedEncryptParty + elapsedEvalCloud + elapsedEvalParty + elapsedDecCloud + elapsedDecParty) / time.Duration(nexperiments)

			elapsedSetupParty_ri = time.Duration(0)
			elapsedSetupParty = time.Duration(0)
			elapsedInputParty = time.Duration(0)
			elapsedEncryptCloud = time.Duration(0)
			elapsedEncryptParty = time.Duration(0)
			elapsedEvalCloud = time.Duration(0)
			elapsedEvalParty = time.Duration(0)
			elapsedDecCloud = time.Duration(0)
			elapsedDecParty = time.Duration(0)
		}

	}

	l.Printf("\t\t\t||                                                                                                                       ||")
	l.Printf("\t\t\t||            - Mean Setup:     %12s                                                                             ||", SetupMean)
	l.Printf("\t\t\t||            - Mean GenInputs: %12s                                                                             ||", GenInputMean)
	l.Printf("\t\t\t||            - Mean EncPhase:  %12s                                                                             ||", EncPhaseMean)
	l.Printf("\t\t\t||            - Mean Agg:       %12s                                                                             ||", AggMean)
	l.Printf("\t\t\t||            - Mean DecPhase:  %12s                                                                             ||", DecMean)
	l.Printf("\t\t\t||            - Mean Total:     %12s                                                                             ||", TotalMean)
	l.Printf("\t\t\t||=======================================================================================================================||")

	return


}
// ----------------------------------------------------------------------------------------------//
// Generates the invidividual secret key "ski" for each Data Owner Party P[i]
// ----------------------------------------------------------------------------------------------//

// genParties initializes each data owner party and generates their individual secret key.
// Inputs:
// - aggring: Aggregation ring structure containing the cryptographic parameters
// - secretkeysampler: Sampler used to generate secret keys
// - NumParties: Total number of data owner parties
// Outputs:
// - Returns an array of *party structures, each containing an initialized secret key for a data owner

func genParties(aggring *AggRings, secretkeysampler ring.Sampler, NumParties int) []*party {

	// Allocate memory for each party's structure and the necessary shares for protocol operations
	P := make([]*party, NumParties)

	// Track the setup time for party initialization
	elapsedSetupParty += runTimedParty(func() {

		// Generate and initialize the secret key for each party
		for i := range P {
			pi := &party{}                     // Create a new party instance
			pi.sk = secretkeysampler.ReadNew() // Generate a new secret key using the provided sampler
			aggring.ringQ.NTT(pi.sk, pi.sk)    // Transform the secret key to the NTT domain for efficient polynomial operations
			P[i] = pi                          // Assign the initialized party to the party array
		}
	}, len(P))

	return P // Return the array of initialized parties
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Generates the individual share "ri" for each Data Owner Party P[i]
// ----------------------------------------------------------------------------------------------//

// genSetupShare initializes the share "ri" for each party so that the sum of all "ri" values is zero.
// Inputs:
// - aggring: Aggregation ring structure containing the cryptographic parameters
// - uniformsampler: Sampler used to generate uniform random values for each "ri"
// - P: Array of parties, each of which will receive an "ri" share
func genSetupShare(aggring *AggRings, uniformsampler ring.Sampler, P []*party) {

	// Create each party and allocate memory for all shares required by the protocol.
	// The goal is to set up the shares "ri" for each party such that:
	//     NTT(MForm(r_{L-1})) = -(NTT(MForm(r_0)) + NTT(MForm(r_1)) + ... + NTT(MForm(r_{L-2})))
	// This ensures that the sum of all shares is zero, which is essential for the aggregation protocol.

	// Track the setup time for generating shares for all parties

	// Initialize the last party's "r" to accumulate the negated sum of previous shares
	P[len(P)-1].r = aggring.ringQ.NewPoly()

	// Generate the "ri" shares for each party except the last
	for _, pi := range P[:len(P)-1] {
		pi.r = uniformsampler.ReadNew() // Generate a new uniform random share for each party
		aggring.ringQ.NTT(pi.r, pi.r)   // Convert the share to the NTT domain
		aggring.ringQ.MForm(pi.r, pi.r) // Apply the Montgomery form transformation

		// Accumulate each share into the last party's share, P[len(P)-1].r, to ensure zero sum
		aggring.ringQ.Add(pi.r, P[len(P)-1].r, P[len(P)-1].r)
	}

	// Negate the accumulated value in the last party's share to achieve a zero-sum across all shares
	aggring.ringQ.Neg(P[len(P)-1].r, P[len(P)-1].r)

}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Generates the input polynomials "mi" for each Data Owner Party P[i] and the expected aggregated result "aggexp"
// ----------------------------------------------------------------------------------------------//

// genInputs initializes the input "mi" for each party and obtains the expected aggregated result aggexp = m1 + m2 + ...
// Inputs:
//   - aggring: Aggregated ring structure containing necessary ring operations and parameters.
//   - lowNormUniformQ: Sampler to generate polynomials with coefficients bounded by a norm.
//   - n: Number of input polynomials to be generated per party.
//   - P: Array of parties participating in the protocol, each having an input "mi" to contribute.
//   - param: Struct containing cryptographic parameters, specifically the modulus level for polynomial generation.
// Outputs:
//   - aggexp: Array of polynomials representing the expected aggregate sum of inputs from all parties.
// ----------------------------------------------------------------------------------------------//

func genInputs(aggring *AggRings, lowNormUniformQ *lowNormSampler, n int, P []*party, param parameters) (aggexp []ring.Poly) {

    elapsedInputParty += runTimedParty(func() {

        aggexp = make([]ring.Poly, n)
        for i := 0; i < n; i++ {
            aggexp[i] = aggring.ringQ.NewPoly() 
        }

        delta := aggring.Delta

        for _, pi := range P {
            pi.input = make([]ring.Poly, n)
            for i := 0; i < n; i++ {
                // Generate in range [-Delta, Delta]
                pi.input[i] = lowNormUniformQ.newPolyLowNorm(delta)

                // Transform to NTT
                aggring.ringQ.NTT(pi.input[i], pi.input[i])

                // Acumulate in aggexp
                aggring.ringQ.Add(pi.input[i], aggexp[i], aggexp[i])
            }
        }

        // Back to standard domain
		for i := 0; i < n; i++ {
            aggring.ringQ.INTT(aggexp[i], aggexp[i]) 
            //aggexp[i].Resize(param.plevel)           
        }

    }, len(P))

    return aggexp
}
// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Generates the input Secret-Key ciphertexts "encInputs" + partial decryptions "partialDec"
// ----------------------------------------------------------------------------------------------//

// encPhase generates encrypted inputs and partial decryptions for each party in the protocol.
// For each input "mi"" of each party, this phase performs an encryption operation on the input
// using the party's secret key ("ski") and random share ("ri"). This process involves sampling a
// random polynomial "a", generating noise "ei", and applying modular transformations to
// secure the data before aggregation. Additionally, each party computes a partial decryption "pdi" share
// that will be used in the final decryption phase.
// Inputs:
//   - aggring: Aggregated ring structure with necessary operations for encryption and decryption.
//   - gaussianSamplerQ: Sampler for generating Gaussian noise polynomials.
//   - uniformSamplerQ: Sampler for generating uniform random polynomials "a" used in encryption.
//   - n: Number of input polynomials to generate per party.
//   - P: Array of parties participating in the protocol, each holding a secret key "sk" and a randomness "r".
//   - param: Struct containing cryptographic parameters, especially the prime modulus level for encryption.
//
// Outputs:
//   - encInputs: Array of encrypted input polynomials for each party, structured as encInputs[party][input].
//   - partialDec: Array of partial decryptions for each party, structured as partialDec[party][input].
func encPhase(aggring *AggRings, gaussianSamplerQ ring.Sampler, uniformSamplerQ ring.Sampler, n int, P []*party) (encInputs [][]ring.Poly, partialDec [][]ring.Poly) { // to point or not to point, that's the question

	// Measure elapsed time for encryption phase across all parties
	elapsedEncryptParty = time.Duration(0)

	tmp := aggring.ringQ.NewPoly() // Temporary polynomial for computation steps

	// Initialize arrays to hold encrypted inputs and partial decryptions for each party
	encInputs = make([][]ring.Poly, len(P))
	partialDec = make([][]ring.Poly, len(P))
	for i := 0; i < len(P); i++ {
		encInputs[i] = make([]ring.Poly, n)
		partialDec[i] = make([]ring.Poly, n)
	}

	FirstNoise := make([][]ring.Poly, len(P))
	SecondNoise := make([][]ring.Poly, len(P))
	for i := 0; i < len(P); i++ {
		FirstNoise[i] = make([]ring.Poly, n)
		SecondNoise[i] = make([]ring.Poly, n)
	}

	// Step 1: Loop over each input index "j" to sample random polynomials "a" for encryption. Sample different "a" per each consecutive correlated encryption (a total of "n")
	for j := 0; j < n; j++ {

		// Encryption: a*(si + ri) + e + Δ*m[i]
		a := uniformSamplerQ.ReadNew() // Sample a uniform random polynomial "a"
		aggring.ringQ.NTT(a, a)        // Convert "a" to NTT domain for efficient operations

		elapsedEncryptParty += runTimedParty(func() {

			for i := 0; i < len(P); i++ { // Loop over each party to encrypt their inputs
				// c = e
				encInputs[i][j] = aggring.ringQ.NewPoly()
				gaussianSamplerQ.AtLevel(aggring.ringQ.MaxLevel()).ReadAndAdd(encInputs[i][j]) // Generate Gaussian noise "e"

				// c = NTT(e)
				aggring.ringQ.NTT(encInputs[i][j], encInputs[i][j]) // Convert noise to NTT domain

				//---------------------Noises---------------------//

				// Copy the noise e to FirstNoise[i][j] for later use
				FirstNoise[i][j] = aggring.ringQ.NewPoly()
				FirstNoise[i][j].Copy(encInputs[i][j])

				SecondNoise[i][j] = aggring.ringQ.NewPoly()
				gaussianSamplerQ.AtLevel(aggring.ringQ.MaxLevel()).ReadAndAdd(SecondNoise[i][j]) // Generate Gaussian noise "e"

				// c = NTT(e)
				aggring.ringQ.NTT(SecondNoise[i][j], SecondNoise[i][j]) // Convert noise e' to NTT domain

				//----------------Noises finish------------------//

				// tmp_1 = NTT(m * Q/P) for scaling. No necessary. The message is already scaled by Delta in inputs.

				// c = NTT(m * (Q/P) + e)
				aggring.ringQ.Add(encInputs[i][j], P[i].input[j], encInputs[i][j]) // Add scaled message to noise

				// NTT(pdi) = NTT(a)*NTT(MF(ski))) = NTT(a*ski)
				partialDec[i][j] = aggring.ringQ.NewPoly()
				aggring.ringQ.MulCoeffsMontgomery(P[i].sk, a, partialDec[i][j])

				//ringQ.MForm(pd[i], pd[i])
				// c = NTT(m * (Q/P) + e) + NTT(a) * MForm(NTT(ri))
				// c = NTT(m * (Q/P) + e + a*ri)
				//ringQ.MForm(c[i], c[i])
				aggring.ringQ.MulCoeffsMontgomery(a, P[i].r, tmp) // Multiply "a" with party's randomness "r" and add to ciphertext
				aggring.ringQ.Add(tmp, encInputs[i][j], encInputs[i][j])

				// c = NTT(m * (Q/P) + e + a*ri + a*ski)
				aggring.ringQ.Add(encInputs[i][j], partialDec[i][j], encInputs[i][j]) // Complete encryption: add a*ski

				//----------------partialdec = NTT(ask + e - e') ----------------------------------------------------//

				aggring.ringQ.Add(partialDec[i][j], FirstNoise[i][j], partialDec[i][j])
				aggring.ringQ.Sub(partialDec[i][j], SecondNoise[i][j], partialDec[i][j])

			}
		}, len(P))
	}

	// Output: Return the encrypted inputs and the partial decryption shares for each party
	return encInputs, partialDec
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Computes the aggregation c0 + c1 + ... + c_{L - 1} = cagg
// ----------------------------------------------------------------------------------------------//

// evalPhase aggregates encrypted inputs from each party into a single encrypted sum.
// This function takes in multiple ciphertexts from each party and sums them to produce
// an aggregate ciphertext `cagg`, representing the total encrypted data of all parties.

// Inputs:
//   - aggring: Aggregated ring structure containing operations required for polynomial arithmetic.
//   - n: The number of encrypted inputs per party.
//   - encInput: A 2D array of encrypted polynomials, where encInput[party][input] represents
//               the encrypted input from each party for each data point.
//   - param: Parameter structure holding encryption settings and the target modulus level.

// Outputs:
//   - encShareAgg: Array of aggregated polynomials, where each entry corresponds to the sum of
//                  encrypted inputs for a given data point across all parties.

func evalPhase(aggring *AggRings, n int, encInput [][]ring.Poly) (encShareAgg []ring.Poly) {

	// Measure and add time spent on this function to "elapsedEvalCloud"
	elapsedEvalCloud += runTimed(func() {

		// Initialize "encShareAgg" to store the aggregated ciphertext for each input
		encShareAgg = make([]ring.Poly, n)
		for i := 0; i < n; i++ {
			encShareAgg[i] = aggring.ringQ.NewPoly() // Allocate space for each aggregate polynomial
		}

		// Aggregation process: iterate over each input index "j"
		for j := 0; j < n; j++ { // Loop through each ciphertext slot
			for i := 0; i < len(encInput); i++ { // Loop through each party
				// Add the encrypted input from party "i" to the aggregate for the current input "j"
				aggring.ringQ.Add(encInput[i][j], encShareAgg[j], encShareAgg[j])
			}

			
		}
	})

	// Return the aggregated ciphertext array
	return encShareAgg
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Gathers all partial decryptions "partialDec" together with "encShareAgg" to obtain the result
// ----------------------------------------------------------------------------------------------//

// decPhase combines all partial decryptions from each party with the aggregated encrypted
// share "encShareAgg" to compute the final decrypted result in parallel.
//
// Inputs:
//   - aggring: Aggregated ring structure containing operations required for polynomial arithmetic.
//   - partialDec: A 2D array of partial decryptions from each party, where partialDec[party][input]
//     represents the partial decryption share for each data point from each party.
//   - encShareAgg: Array of aggregated encrypted polynomials, representing the sum of encrypted inputs
//     across all parties for a given data point.
//   - n: The number of encrypted inputs per party.
//   - param: Parameter structure holding encryption settings and the target modulus level.
//
// Outputs:
//   - recShare: Array of decrypted polynomials, where each entry corresponds to the final decrypted
//     result after combining all partial decryptions and the aggregated encrypted share.

// ----------------------------------------------------------------------------------------------//
// Gathers all partial decryptions "partialDec" together with "encShareAgg" to obtain the result
// ----------------------------------------------------------------------------------------------//
func decPhase(aggring *AggRings, partialDec [][]ring.Poly, encShareAgg []ring.Poly, n int, param parameters, P []*party) (recShare []ring.Poly) {

	// Assumes that messages are masked with known randomness by all input parties,
	// allowing the cloud to securely perform decryption.

	// Initialize a timer to measure the decryption phase duration.
	elapsedDecCloud = time.Duration(0)
	elapsedDecCloud += runTimedParty(func() {

		//buff := aggring.ringQ.NewPoly()

		// Preparation: Initialize data structures
		// Polynomials to store the decrypted results.
		recShare = make([]ring.Poly, n)

		for j := 0; j < n; j++ { // running through ciphertexts per party

			recShare[j] = aggring.ringQ.NewPoly()

			// Aggregate shares
			// Store in recShare the content of encShareAgg
			aggring.ringQ.Add(encShareAgg[j], recShare[j], recShare[j])
			for i := 0; i < len(partialDec); i++ {
				aggring.ringQ.Sub(recShare[j], partialDec[i][j], recShare[j])
			}

			// Eliminate NTT 
			aggring.ringQ.INTT(recShare[j], recShare[j])


		}

	}, len(P))

	// Return the decrypted results.
	return recShare
}









// DivideByBigInt divides each element of a float64 slice by the value of a big.Int divisor.
func DivideByBigInt(slice []float64, divisor *big.Int) []float64 {
	result := make([]float64, len(slice))
	divisorFloat, _ := new(big.Float).SetInt(divisor).Float64()
	for i, v := range slice {
		result[i] = v / divisorFloat
	}
	return result
}

// ----------------------------------------------------------------------------------------------//

func PolyToSignedFloat64Slice(pol ring.Poly, r *ring.Ring) []float64 {
	level := pol.Level()
	N := r.N()

	// 1. Calcular el módulo Q total hasta el nivel actual del polinomio
	Q := big.NewInt(1)
	for _, qi := range r.ModuliChain()[:level+1] {
		Q.Mul(Q, new(big.Int).SetUint64(qi))
	}

	// 2. Calcular Q/2 para el centrado
	QHalf := new(big.Int).Div(Q, big.NewInt(2))

	bigCoeffs := make([]*big.Int, N)
	for i := range bigCoeffs {
		bigCoeffs[i] = new(big.Int)
	}

	// 3. Reconstruir el BigInt a partir del sistema RNS
	r.AtLevel(level).PolyToBigint(pol, 1, bigCoeffs)

	coeffs := make([]float64, N)
	for i := 0; i < N; i++ {
		// 4. Centrado: Si el valor es > Q/2, es un número negativo
		if bigCoeffs[i].Cmp(QHalf) > 0 {
			bigCoeffs[i].Sub(bigCoeffs[i], Q)
		}

		f, _ := new(big.Float).SetInt(bigCoeffs[i]).Float64()
		coeffs[i] = f
	}
	return coeffs
}
