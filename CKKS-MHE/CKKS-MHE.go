// Copyright 2026 Miguel Morona-Mínguez, Alberto Pedrouzo-Ulloa
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
// limitations under the License.

package main

import (
	"log"
	"math"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// ----------------------------------------------------------------------------------------------//
// For more details see the report discussing the implementation runtimes of aggregation with baseline BFV and the BFVbaselineAgg folder in "https://github.com/apedrouzoulloa/mkagg"
// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// CKKS-MHE.go can be run as: go run ./CKKS-MHE.go NumParties Goroutines NumModelParameters
// - By default: NumParties = 32, Goroutines = 1, NumModelParameters = 8192000
//
// Cryptosystem parameters are defined in "paramsDef"
// - By default:
//	paramsDef := ckks.ParametersLiteral{ //CKKSRealParamsN13QP218
//		LogN: 13,
//		LogQ: []int{60, 60}, // Q ~ 120 bits
//		LogDefaultScale: 118, // 2^118
//	}
//
// ----------------------------------------------------------------------------------------------//
// Comments: The code of this script was started by relying on the Examples folder available in Lattigo "https://github.com/tuneinsight/lattigo"
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

// To measure time execution of the same operations which are simultaneously done in "parallel" by L parties
func runTimedParty(f func(), L int) time.Duration {
	start := time.Now()
	f()
	return time.Duration(time.Since(start).Nanoseconds() / int64(L))
}

// party: defines the inputs and private information of party P_i
type party struct {
	sk                       *rlwe.SecretKey              // sk_i (individual secret key of party P_i)
	ckgShare                 multiparty.PublicKeyGenShare // Party share used during the generation of the collective public key
	input                    []float64                    // Unidimensional array containing the vector associated with the model update
	PaddedNumModelParameters int                          // Size of the vector of model updates with zero padding (till the next multiple of params.N())
}

// multTask: structure used for the parallelization with Go routines in aggregation phase
type multTask struct {
	wg              *sync.WaitGroup
	op1             []*rlwe.Ciphertext
	op2             []*rlwe.Ciphertext
	res             []*rlwe.Ciphertext
	elapsedmultTask time.Duration
}

// Definition of variables used to measure runtime for specific steps of the aggregation protocol
var elapsedEncryptParty time.Duration
var elapsedEncryptCloud time.Duration
var elapsedSetupParty time.Duration
var elapsedCKGCloud time.Duration
var elapsedCKGParty time.Duration
var elapsedPCKSCloud time.Duration
var elapsedPCKSParty time.Duration
var elapsedEvalCloudCPU time.Duration
var elapsedEvalCloud time.Duration
var elapsedEvalParty time.Duration
var elapsedDecCloud time.Duration
var elapsedDecParty time.Duration
var elapsedHomKeySwitchCloud time.Duration
var elapsedHomKeySwitchParty time.Duration

var elapsedDecodeParty time.Duration
var elapsedDecodeCloud time.Duration

var SetupMean time.Duration
var GenInputMean time.Duration
var EncPhaseMean time.Duration
var AggMean time.Duration
var DecMean time.Duration
var TotalMean time.Duration

// main function. Example executions:
// go run ./CKKS-MHE.go 32 1 8192000
func main() {

	var err error

	file, err := os.Create("ckks-results-flta.txt")
	check(err)
	defer file.Close() 
	l := log.New(file, "", 0)

	
	//l := log.New(os.Stderr, "", 0)

	// Default number of input parties. NumParties
	L := 32
	

	// If at least one argument is provided, attempt to parse the first argument as an integer
	// to redefine the number of parties (L).
	if len(os.Args[1:]) >= 1 {
		L, err = strconv.Atoi(os.Args[1])
		check(err)
	}

	// Default number of goroutines to use.
	NGoRoutine := 1

	// If a second argument is provided, attempt to parse it as an integer
	// to redefine the number of goroutines.
	if len(os.Args[1:]) >= 2 {
		NGoRoutine, err = strconv.Atoi(os.Args[2])
		check(err)
	}

	// Default number of model parameters
	NumModelParameters := 1000 //8192000 // This equals 2^13 * 1000 (2000ctx as maxslots is 2^12) (As there is only aggregation, we can suppone a conjugate invariant ring) //1024 total

	// If a third argument is provided, attempt to parse it as an integer
	// to redefine the number of model parameters.
	if len(os.Args[1:]) >= 3 {
		NumModelParameters, err = strconv.Atoi(os.Args[3])
		check(err)
	}

	paramsDef := ckks.ParametersLiteral{ //CKKSRealParamsN13QP218
		LogN: 13,
		LogQ: []int{60, 60}, // Q ~ 120 bits
		//LogP:            []int{35},
		LogDefaultScale: 118, // 2^118
		//RingType:        ring.ConjugateInvariant,
	}

	nexperiments := 25

	for z := 0; z < nexperiments; z++ {

		l.Printf("Number of experiment: n=%2d", z)

		// Si Log(B) ~ 100 bits, la precisión es de ~20 bits y el q es de 121 bits

		// Initialize cryptographic parameters from "paramsDef"
		params, err := ckks.NewParametersFromLiteral(paramsDef)
		check(err)

		// Create a cryptographic random source (CRS) with a keyed PRNG using a predefined seed.
		crs, err := sampling.NewKeyedPRNG([]byte{'b', 'a', 's', 'e', 'l', 'i', 'n', 'e'})
		check(err)

		// Initialize the encoder for encoding and decoding plaintexts
		encoder := ckks.NewEncoder(params)

		// Generate target secret and public key pair for testing decryption and correctness
		tsk, tpk := rlwe.NewKeyGenerator(params).GenKeyPairNew()

		// Create each party and allocate the memory for all the shares needed in the protocol
		P := genparties(params, L)
		l.Println("> Initialization of Parties")

		// Start the process, only 1 aggregation round is executed

		// Generate inputs and calculate the expected result for verification purposes.
		// The inputs are padded with zeros for proper encoding.

		expRes := genInputs(params, P, NumModelParameters, 1) // BoundInputs = 1

		l.Printf("> Input generation\n \tNum parties: %d, PaddedwithZeros-NumModelParameters: %d\n", len(P), len(expRes))

		// 1) Collective public key generation: Perform collective key generation (CKG) to compute a public key shared by all parties.
		pk := ckgphase(params, crs, P, l)

		l.Printf("\tSetup done (cloud: %s, party: %s). Total: %s\n", elapsedCKGCloud, elapsedCKGParty+elapsedSetupParty, elapsedCKGCloud+elapsedCKGParty+elapsedSetupParty)

		// 2) Encryption phase: Each party encrypts their inputs using the shared public key.
		encInputs := encPhase(params, P, pk, encoder, l)

		// 3) Evaluation phase: Perform homomorphic aggregation on the encrypted inputs.
		encRes := evalPhase(params, NGoRoutine, encInputs, l)
		encInputs = nil

		// 4) Public Collective Key Switching (PCKS): Parties collaboratively switch the encrypted result to the target secret key.
		encOut := pcksPhase(params, tpk, encRes, P, l)
		encRes = nil
		//P = nil
		l.Printf("Size of result\t: Number of ciphertexts: %d ciphertexts\n", len(encOut))

		// Decryption phase: Use the target secret key (tsk) to decrypt the final result.
		l.Println("> Decrypt Phase")

		// We do the decryption step diectly in the code, avoiding the use of an external function

		decryptor := rlwe.NewDecryptor(params, tsk)

		// Prepare plaintexts to hold decrypted data.
		ptres := make([]*rlwe.Plaintext, len(encOut))
		for i := range encOut {
			ptres[i] = ckks.NewPlaintext(params, params.MaxLevel())
		}

		// Only 1 decryption is run, but would be done by all parties who know the tsk. So we use runTimed instead of runTimedParty to measure the runtime
		elapsedDecParty = runTimed(func() {
			for i := range encOut {
				decryptor.Decrypt(encOut[i], ptres[i]) // Decrypt each ciphertext into plaintext
				encOut[i] = nil
			}
		})

		elapsedDecCloud = time.Duration(0)
		l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedDecCloud, elapsedDecParty)
		
		// Decode phase: convert plaintext type into float.
		l.Println("> Result and decode phase:")

		// Check the result
		// Decode plaintexts to retrieve the final result as a list of floats.
		res := make([]float64, len(expRes))

		partialRes := make([]float64, params.MaxSlots()) // Temporary buffer for decoding.
		
		elapsedDecodeParty = runTimed(func() {
		for i := range ptres {
			check(encoder.Decode(ptres[i], partialRes)) // Decode plaintext to float64 values.
			ptres[i] = nil
			for j := range partialRes {
				res[(i*len(partialRes) + j)] = partialRes[j] // Copy values to the final result array.
			}
		}})
		
		elapsedDecodeCloud = time.Duration(0)
		l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedDecodeCloud, elapsedDecodeParty)		

		MaxDiff := 0.0

		// Validate the result by comparing it against the expected result.
		l.Printf("\t%v\n", res[:10])
		l.Printf("\t%v\n", expRes[:10])
		for i := range expRes {
			if math.Abs(expRes[i]-res[i]) > MaxDiff {
				MaxDiff = math.Abs(expRes[i] - res[i])
			}
			if math.Abs(expRes[i]-res[i]) > math.Exp2(float64(-20)) { //Epsilon = 2^-20
				// Log error details if there is a mismatch.
				l.Printf("\tincorrect\n first error in position [%d]\n", i)
				l.Printf("Max Error: %v", MaxDiff)
				l.Printf("Epsilon: %v", math.Exp2(float64(-20)))
				l.Printf("> Finished (total cloud: %s, total party: %s)\n", elapsedCKGCloud+elapsedEncryptCloud+elapsedEvalCloud+elapsedPCKSCloud, elapsedCKGParty+elapsedEncryptParty+elapsedEvalParty+elapsedPCKSParty+elapsedDecParty)
				return
			}
		}
		// If all results match, log success.
		l.Println("\tcorrect")
		l.Printf("Max Error: %v", MaxDiff)
		l.Printf("Epsilon: %v", math.Exp2(float64(-20)))
		l.Printf("> Finished (total cloud: %s, total party: %s)\n", elapsedCKGCloud+elapsedEncryptCloud+elapsedEvalCloud+elapsedPCKSCloud+elapsedDecCloud, elapsedCKGParty+elapsedEncryptParty+elapsedEvalParty+elapsedPCKSParty+elapsedDecParty)

		SetupMean += (elapsedCKGCloud + elapsedCKGParty + elapsedSetupParty) / time.Duration(nexperiments)
		EncPhaseMean += (elapsedEncryptCloud + elapsedEncryptParty) / time.Duration(nexperiments)
		AggMean += (elapsedEvalCloud + elapsedEvalParty) / time.Duration(nexperiments)
		DecMean += (elapsedPCKSCloud + elapsedPCKSParty + elapsedDecCloud + elapsedDecParty + elapsedDecodeCloud + elapsedDecodeParty) / time.Duration(nexperiments)
		TotalMean += (elapsedCKGCloud + elapsedCKGParty + elapsedSetupParty + elapsedEncryptCloud + elapsedEncryptParty + elapsedEvalCloud + elapsedEvalParty + elapsedPCKSCloud + elapsedPCKSParty + elapsedDecCloud + elapsedDecParty) / time.Duration(nexperiments)

		elapsedSetupParty = time.Duration(0)
		elapsedEncryptCloud = time.Duration(0)
		elapsedEncryptParty = time.Duration(0)
		elapsedEvalCloud = time.Duration(0)
		elapsedEvalParty = time.Duration(0)
		elapsedDecCloud = time.Duration(0)
		elapsedDecParty = time.Duration(0)
		elapsedPCKSCloud = time.Duration(0)
		elapsedPCKSParty = time.Duration(0)
		elapsedDecodeCloud = time.Duration(0)
		elapsedDecodeParty = time.Duration(0)		

	}

	l.Printf("\t\t\t||                                                                                                                       ||")
	l.Printf("\t\t\t||=======================================================================================================================||")
	l.Printf("\t\t\t||            - Mean Setup:     %12s                                                                             ||", SetupMean)
	l.Printf("\t\t\t||            - Mean EncPhase:  %12s                                                                             ||", EncPhaseMean)
	l.Printf("\t\t\t||            - Mean Agg:       %12s                                                                             ||", AggMean)
	l.Printf("\t\t\t||            - Mean DecPhase:  %12s                                                                             ||", DecMean)
	l.Printf("\t\t\t||            - Mean Total:     %12s                                                                             ||", TotalMean)

}

// ----------------------------------------------------------------------------------------------//
// Generates the individual secret key and input model updates of size "NumModelParameters"
// for each Data Owner Party P[i]
// ----------------------------------------------------------------------------------------------//

// genparties initializes each data owner party and generates their individual secret key.
// Inputs:
// - params: Cryptographic parameters necessary for key generation
// - N: Total number of data owner parties
// Outputs:
// - Returns an array of *party structures, each containing an initialized secret key for a data owner

func genparties(params ckks.Parameters, NumParties int) []*party {

	// Allocate memory for each party's structure and the necessary shares for protocol operations
	P := make([]*party, NumParties)

	// Track the setup time for party initialization
	elapsedSetupParty += runTimedParty(func() {

		// Initialize each party and generate their individual secret keys
		for i := range P {
			pi := &party{}                                         // Create a new party instance
			pi.sk = rlwe.NewKeyGenerator(params).GenSecretKeyNew() // Generate a new secret key using the provided parameters
			P[i] = pi                                              // Assign the initialized party to the party array
		}
	}, NumParties)

	return P // Return the array of initialized parties
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Generates the inputs for each data owner party P[i] and calculates the expected result.
// The inputs are randomized values constrained by the provided BoundInputs and the parameters of the model.
// ----------------------------------------------------------------------------------------------//

// genInputs initializes the input values for each party and computes the expected results based on
// the inputs from all parties.
// Inputs:
// - params: Cryptographic parameters required for padding and modular operations
// - P: Array of data owner parties, each containing an input for the model
// - NumModelParameters: Total number of model parameters (used to determine input padding size)
// - BoundInputs: Upper bound for random input generation
// Outputs:
// - expRes: A slice of uint64 values representing the aggregated expected result of all parties' inputs

func genInputs(params ckks.Parameters, P []*party, NumModelParameters int, BoundInputs float64) (expRes []float64) {

	// Generate Inputs for each party
	for _, pi := range P {
		// Determine the number of model parameters to pad based on polynomial degree and model parameters
		if params.MaxSlots() >= NumModelParameters { // If polynomial degree is greater than or equal to the number of model parameters, no padding is needed
			pi.PaddedNumModelParameters = params.MaxSlots()
		} else { // If polynomial degree is less than the number of model parameters, calculate required padding
			pi.PaddedNumModelParameters = int(math.Ceil(float64(NumModelParameters)/float64(params.MaxSlots()))) * params.MaxSlots()
		}

		pi.input = make([]float64, pi.PaddedNumModelParameters)

		for i := range pi.input {
			if i < NumModelParameters {
				pi.input[i] = (sampling.RandFloat64(-BoundInputs/float64(len(P)), BoundInputs/float64(len(P)))) //(min/N,max/N)

			} else {
				pi.input[i] = 0
			}
		}
	}

	// Allocate memory for the expected result array
	expRes = make([]float64, P[0].PaddedNumModelParameters)

	// Generate the Aggregation Expected Results by summing inputs from all parties

	for _, pi := range P {
		for i := range pi.input {
			expRes[i] += pi.input[i]
		}
	}

	// Return the aggregated expected result
	return
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Conducts the Collective Key Generation (CKG) phase where each party generates a share of the
// public key and aggregates the shares to form the final public key.
// ----------------------------------------------------------------------------------------------//

// ckgphase performs the Collective Key Generation (CKG) protocol, where each party generates
// its share of the public key and the shares are then combined to create the final public key.
// Inputs:
// - params: Cryptographic parameters required for key generation
// - crs: It is used for generating the common random polynomial
// - P: Array of data owner parties, each generating a share of the public key
// Outputs:
// - Returns the aggregated public key generated from all party shares

func ckgphase(params ckks.Parameters, crs sampling.PRNG, P []*party, l *log.Logger) *rlwe.PublicKey {

	//l := log.New(os.Stderr, "", 0)

	// Log the start of the CKG phase
	l.Println("> CKG Phase")

	// Initialize the public key generation protocol
	ckg := multiparty.NewPublicKeyGenProtocol(params)

	// Allocate a combined share for the public key, initially empty
	ckgCombined := ckg.AllocateShare()

	// Allocate space for each party's share in the key generation
	for _, pi := range P {
		pi.ckgShare = ckg.AllocateShare()
	}

	var crp multiparty.PublicKeyGenCRP // added by APU

	// Record the time taken for the party-side share generation phase
	elapsedCKGParty = runTimedParty(func() {

		// Sample the crp from the crs
		crp = ckg.SampleCRP(crs)

		// Each party generates its own share of the public key based on its secret key and the CRP
		for _, pi := range P {
			ckg.GenShare(pi.sk, crp, &pi.ckgShare)
		}
	}, len(P))

	// Create a new public key to hold the final result
	pk := rlwe.NewPublicKey(params)

	// Record the time taken for the cloud-side aggregation and public key generation
	elapsedCKGCloud = runTimed(func() {

		// Aggregate each party's share into the combined share
		for _, pi := range P {
			ckg.AggregateShares(pi.ckgShare, ckgCombined, &ckgCombined)
		}

		// Generate the final public key by using the combined share and the CRP
		ckg.GenPublicKey(ckgCombined, crp, pk)
	})

	// Log the time spent on the cloud and party-side operations
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedCKGCloud, elapsedCKGParty)

	// Return the final generated public key
	return pk
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Executes the encryption phase where each party encrypts its input data into ciphertexts using
// the public key. Each party’s input is divided into multiple ciphertexts based on the model's
// parameter size, and encryption is performed using the provided encoder.
// ----------------------------------------------------------------------------------------------//

// encPhase performs the encryption of each party's input data, where each input is split into
// multiple ciphertexts depending on the model's parameter size and the party's input vector.
// Inputs:
// - params: Cryptographic parameters necessary for the encryption
// - P: Array of data owner parties, each containing input data to be encrypted
// - pk: Public key used for encrypting the input data
// - encoder: Encoder used for encoding the input data into plaintexts before encryption
// Outputs:
// - encInputs: A 2D slice of ciphertexts, where each party's encrypted inputs are stored in ciphertexts

func encPhase(params ckks.Parameters, P []*party, pk *rlwe.PublicKey, encoder *ckks.Encoder, l *log.Logger) (encInputs [][]*rlwe.Ciphertext) {

	//l := log.New(os.Stderr, "", 0)

	// Determine the number of ciphertexts each party needs to generate based on the number of model parameters
	NumCiphertextsPerParty := int(math.Ceil(float64(P[0].PaddedNumModelParameters) / float64(params.MaxSlots())))

	//inputs := CombineInputs(P)
	//blocksInputs := ChunkAndPadInputsFloat64(inputs, params)

	// Initialize the encInputs 2D array to hold ciphertexts for each party and each model parameter
	// encInputs[i][j], i through Parties, j through Model parameters
	encInputs = make([][]*rlwe.Ciphertext, len(P))
	for i := range encInputs {
		encInputs[i] = make([]*rlwe.Ciphertext, NumCiphertextsPerParty)
	}

	for i := range encInputs {
		for k := range encInputs[i] {
			encInputs[i][k] = ckks.NewCiphertext(params, 1, params.MaxLevel()) // Allocate a new ciphertext for each entry
		}
	}

	// Start the encryption phase
	l.Println("> Encrypt Phase")

	// Create an encryptor using the public key
	encryptor := rlwe.NewEncryptor(params, pk)

	// Create a plaintext object and an array to hold the input values for encryption
	pt := ckks.NewPlaintext(params, params.MaxLevel())
	elapsedEncryptParty = runTimedParty(func() {
		for i, pi := range P {
			for k := 0; k < NumCiphertextsPerParty; k++ {
				//encoder.Encode(EmbedFloat64SliceToComplex(pi.input[(k*params.N()):((k+1)*params.N())]), pt) // es uno menos en la segunda parte porque go indexa [0:n] como los valores 0, 2, ..., n - 1
				check(encoder.Encode(pi.input[(k*params.MaxSlots()):((k+1)*params.MaxSlots())], pt)) // Encode the input into the plaintext
				check(encryptor.Encrypt(pt, encInputs[i][k]))                                        // Encrypt the plaintext and store the ciphertext
			}

		}
	}, len(P))

	elapsedEncryptCloud = time.Duration(0)
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedEncryptCloud, elapsedEncryptParty)

	// Return the array of ciphertexts generated by all parties
	return
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Executes the evaluation phase where encrypted model updates are combined across multiple layers.
// This phase performs additions of encrypted ciphertexts in parallel using multiple Go routines.
// ----------------------------------------------------------------------------------------------//

// evalPhase performs the evaluation of encrypted model updates in multiple layers. In each layer,
// ciphertexts from different parties are added together in parallel, and the result is propagated
// to the next layer. The task is split among multiple Go routines to improve performance.
// Inputs:
// - params: Cryptographic parameters necessary for the evaluation
// - NGoRoutine: Number of Go routines to use for parallel computation
// - encInputs: A 2D slice of ciphertexts, where each party's encrypted inputs are stored
// Outputs:
// - encRes: A slice of ciphertexts representing the final evaluation result

func evalPhase(params ckks.Parameters, NGoRoutine int, encInputs [][]*rlwe.Ciphertext, l *log.Logger) (encRes []*rlwe.Ciphertext) {

	//l := log.New(os.Stderr, "", 0)

	// Determine the number of layers needed for the evaluation, based on the number of inputs
	SizeEncLayers := make([]int, 0)
	SizeEncLayers = append(SizeEncLayers, len(encInputs))
	endingLayersStructure := 0 // // endingLayerStructure is set to 1 when only 1 layer remains (i.e., nLayer = 1)
	for nLayer := (len(encInputs)/2 + (len(encInputs) & 1)); (nLayer > 0) && (endingLayersStructure == 0); nLayer = ((nLayer >> 1) + (nLayer & 1)) {
		if nLayer == 1 {
			endingLayersStructure = 1
		}
		SizeEncLayers = append(SizeEncLayers, nLayer)
	}

	// Create a channel for task distribution among Go routines and a WaitGroup for synchronization
	tasks := make(chan *multTask)
	workers := &sync.WaitGroup{}
	workers.Add(NGoRoutine)

	// Launch Go routines to process tasks in parallel
	for i := 1; i <= NGoRoutine; i++ {
		go func(i int) {

			// Un evaluador independiente para esta goroutine
			evaluator := ckks.NewEvaluator(params, nil)

			for task := range tasks {
				task.elapsedmultTask = runTimed(func() {
					for indCiphertext := range task.op1 {
						evaluator.Add(
							task.op1[indCiphertext],
							task.op2[indCiphertext],
							task.res[indCiphertext],
						)
					}
				})
				task.wg.Done()
			}
			workers.Done()
		}(i)
	}

	// Start the evaluation tasks
	taskList := make([]*multTask, 0)
	l.Println("> Eval Phase")

	// Scale and shift values used to divide the input layers for processing
	scale := 2
	shift := 1

	// Execute the evaluation phase by adding encrypted ciphertexts across layers
	elapsedEvalCloud = runTimed(func() {
		for i, layer := range SizeEncLayers[:len(SizeEncLayers)-1] {

			nextLayer := SizeEncLayers[i+1]
			l.Println("\tEncrypted model updates added in layer", i, ":", layer, "->", nextLayer)
			wg := &sync.WaitGroup{}

			wg.Add(layer / 2) // Each layer will add pairs of ciphertexts

			for j := 0; j < nextLayer; j++ {
				// Skip certain tasks based on the layer size and index
				if !((2 * nextLayer) > layer) || !(j == (nextLayer - 1)) {
					// Create a task to add the ciphertexts from two input layers
					task := multTask{wg, encInputs[scale*j], encInputs[scale*j+shift], encInputs[scale*j], 0}
					taskList = append(taskList, &task)
					tasks <- &task // Send the task to be processed by a Go routine
				}
			}
			wg.Wait() // Wait for all tasks in the current layer to finish
			scale = 2 * scale
			shift = 2 * shift
		}
	})

	// Compute total time taken for the aggregator-side processing of the tasks
	elapsedEvalCloudCPU = time.Duration(0)
	for _, t := range taskList {
		elapsedEvalCloudCPU += t.elapsedmultTask
	}

	// There is no computation done by parties in this phase, so elapsedEvalParty is 0
	elapsedEvalParty = time.Duration(0)
	l.Printf("\tdone (cloud: %s (wall: %s), party: %s)\n",
		elapsedEvalCloudCPU, elapsedEvalCloud, elapsedEvalParty)

	// Close the task channel and wait for all workers to complete
	close(tasks)
	workers.Wait()

	// The final result of the evaluation is stored in the first element of the input array
	encRes = encInputs[0]
	encInputs = nil

	// Return the array of output ciphertexts generated by the aggregator
	return
}

// ----------------------------------------------------------------------------------------------//

// ----------------------------------------------------------------------------------------------//
// Performs the Public Collective Key Switching (PCKS) phase, where the global secret key
// is switched to the target public key for decryption. This phase updates encrypted results
// with new keys for decryption using the target public key.
// ----------------------------------------------------------------------------------------------//

// pcksPhase executes the collective public-key switching protocol (PCKS) to transform encrypted results
// under a collective secret key into ciphertexts that can be decrypted with the target secret key.
// The protocol is run in two phases: first, each party computes a share of the key switching operation,
// and then the cloud aggregates the shares and performs the final key switching operation.
// Inputs:
// - params: Cryptographic parameters required for the key switching operation
// - tpk: Target public key, whose corresponding secret key will be used for decryption
// - encRes: A slice of encrypted results that need to be switched to the target key
// - P: A slice of parties that hold the secret keys and will participate in the key switching
// Outputs:
// - encOut: A slice of ciphertexts after the key switching has been applied

func pcksPhase(params ckks.Parameters, tpk *rlwe.PublicKey, encRes []*rlwe.Ciphertext, P []*party, l *log.Logger) (encOut []*rlwe.Ciphertext) {

	//l := log.New(os.Stderr, "", 0)

	// To reduce the use of memory: only two pcksShare and pcksCombined components are used. In practice, both should be an array of size NumParties and pcksShare should change for each party

	// Log that the PCKS phase is starting
	l.Println("> PCKS Phase")

	// Initialize the public key switch protocol with a discrete Gaussian distribution for sampling
	pcks, err := multiparty.NewPublicKeySwitchProtocol(params, ring.DiscreteGaussian{Sigma: 1 << 64, Bound: 6 * (1 << 64)})
	//pcks, err := multiparty.NewPublicKeySwitchProtocol(params, ring.DiscreteGaussian{Sigma: 3.19, Bound: 20})
	check(err)

	// Allocate shares for the key switching operation (using a combined share for aggregation)
	pcksCombined := pcks.AllocateShare(params.MaxLevel()) // emulated protocol
	pcksShare := pcks.AllocateShare(params.MaxLevel())    // emulated protocol

	// Loop over each encrypted result (encRes) that needs to be switched to the target public key
	for i := range encRes {

		// Reallocate a new combined share for each encrypted result
		pcksCombined = pcks.AllocateShare(params.MaxLevel())

		// For each party, generate its share of the key switching operation
		for _, pi := range P {

			// Generate key switching share from the party's secret key
			elapsedPCKSParty += runTimedParty(func() {
				// Generate the share using the party's secret key and the target public key
				pcks.GenShare(pi.sk, tpk, encRes[i], &pcksShare) // "emulated protocol"
			}, len(P))

			// Aggregate the key switching shares in the aggregator side
			elapsedPCKSCloud += runTimed(func() {
				pcks.AggregateShares(pcksShare, pcksCombined, &pcksCombined) // "emulated protocol"

			})
		}

		// Perform the key switching operation to update the ciphertext with the new key
		elapsedPCKSCloud += runTimed(func() {
			pcks.KeySwitch(encRes[i], pcksCombined, encRes[i]) // Perform key switching on the ciphertext
		})

	}

	// The output of the phase is the modified ciphertexts (encRes) after the key switching
	encOut = encRes

	// Log the time taken for the cloud and party operations in the key switching phase
	l.Printf("\tdone (cloud: %s, party: %s)\n", elapsedPCKSCloud, elapsedPCKSParty)

	// Return the array of updated ciphertexts
	return

}

//----------------------------------------------------------------------------------------------//

func EmbedFloat64SliceToComplex(input []float64) []complex128 {
	complexSlice := make([]complex128, len(input))
	for i := 0; i < len(input); i++ {
		complexSlice[i] = complex(input[i], 0)
	}
	return complexSlice
}
