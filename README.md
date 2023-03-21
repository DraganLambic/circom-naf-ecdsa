# Batch Signature Verification

Implementation of batch signature verification in circom. We build this project to prove that you know valid signatures for n messages and n corresponding public keys. ECDSA and ED25519 will be suppored in this work.

## ECDSA
We use a lot of the utilities created by [circom-ecdsa](https://github.com/0xparc/circom-ecdsa).

## TODO
We need to transfer the public keys to blochchain address.

A new idea about ff should to be impled. 

In [circom-ecdsa](https://github.com/0xparc/circom-ecdsa), to do a 256 bit secp256k1 scalar mult, we need 256 times ec_double (2946 constraints) and ec_add (2515 constraints), so it takes 1.3m constraints for a scalar mult.

NAF will be use in the new secp256k1 scalar mult, so it will reduce the constraints from 1.3m to 0.95m for a scalar mult.

In NAF, windows is set to be 4. So we get scalar=sum(a_i*(2**4)^i). i is 0 to 63. so we need 256 tims ec_double and 64+15 ec_add.

## Preparatory Work
We rely on snarkjs, circom to get the job done. And we test our work on AWS r5.8xlarge instance with 32-core 3.1GHz, 256G RAM machine with 1T hard drive and 400G swap. 

C++ witness generator will not work on computes with Apple Silicon (e.g., M1) chips due to some assembly incompatibility.

Please follow the steps [Best Practices for Large Circuits](https://hackmd.io/V-7Aal05Tiy-ozmzTGBYPA?view#Setup-from-scratch). Record NODE_PATH, SNARKJS_PATH, RAPIDSNARK_PATH.

You'll need to download a Powers of Tau file with 2^25 constraints and copy it to data directory. You can get in from [this repository](https://github.com/iden3/snarkjs#7-prepare-phase-2).

## How To Test
Modify NODE_PATH, SNARKJS_PATH, RAPIDSNARK_PATH in build_test.sh and run it! 

it will verify 2 ecdsa signature in zk.
