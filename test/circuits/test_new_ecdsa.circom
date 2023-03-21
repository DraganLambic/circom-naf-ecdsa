pragma circom 2.0.3;

include "../../circuits/new_ecdsa.circom";

component main {public [r,s,msghash,pubkey]} = New_ECDSAVerifyNoPubkeyCheck(64, 4);
//component main {public [in]} = ChosenAddPoint();