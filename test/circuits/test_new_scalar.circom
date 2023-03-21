pragma circom 2.0.3;

include "../../circuits/new_secp256k1.circom";

component main {public [scalar, point]} = Secp256k1ScalarNewMult(64, 4);
//component main {public [in]} = ChosenAddPoint();