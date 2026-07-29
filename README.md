# Pushing-Forward-Multi-Secret-Key-Homomorphic-Encryption-for-Private-Average-Aggregation
GOLANG code of the implementations runtimes on the paper. The codes implement secure aggregation using two state of the art HE schemes: **BFV** & **CKKS**. The names of the folders are the following:

* **CKKS-MHE**
* **BFV-MHE**
* **CKKS-PROPOSED**
* **BFV-PROPOSED**

Each folder is an independent Go module (`go.mod`, `go.sum`, `main.go`).

We recommend to check the [Lattigo](https://github.com/tuneinsight/lattigo) Library, where we compare with MHE and we use the cryptographic primitives. Our `PROPOSED` codes run Lattigo v.5, which is not the latest version but works fine. On the other hand, `MHE` codes run the version v.6.

## Citation
xxx

## Contact
To contact us, you can send an email to mmorona@gts.uvigo.es, apedrouzo@gts.uvigo.es or fperez@gts.uvigo.es
