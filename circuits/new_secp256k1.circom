pragma circom 2.1.6;

include "circom-ecdsa/circuits/ecdsa.circom";
include "circom-ecdsa/circuits/secp256k1.circom";


template Secp256k1Double1(n, k) {
    assert(n == 64 && k == 4);

    signal input in[2][k];
    signal input inbackup[2][k]; //used if in==0

    signal output out[2][k];

    signal tmp[2][k];

    component isXZero = BigIsEqual(4);
    for(var i = 0; i < k; i++){
        isXZero.in[0][i] <== in[0][i];
        isXZero.in[1][i] <== 0;
    }


    component isYZero = BigIsEqual(4);
    for(var i = 0; i < k; i++){
        isYZero.in[0][i] <== in[1][i];
        isYZero.in[1][i] <== 0;
    }
 

    signal NotZero <== 1 - isXZero.out*isYZero.out;  //either 0 or 1

    var x1[4];
    var y1[4];
    for(var i=0;i<4;i++){
        x1[i] = (in[0][i]-inbackup[0][i])*NotZero + inbackup[0][i];
        y1[i] = (in[1][i]-inbackup[1][i])*NotZero + inbackup[1][i];
    }

    var tmp1[2][100] = secp256k1_double_func(n, k, x1, y1);
    for(var i = 0; i < k;i++){
        tmp[0][i] <-- tmp1[0][i];
        tmp[1][i] <-- tmp1[1][i];
    }

    component point_on_tangent = Secp256k1PointOnTangent();
    for(var i = 0; i < k; i++){
        point_on_tangent.x1[i] <== x1[i];
        point_on_tangent.y1[i] <== y1[i];
        point_on_tangent.x3[i] <== tmp[0][i];
        point_on_tangent.y3[i] <== tmp[1][i];
    }
 
    component point_on_curve = Secp256k1PointOnCurve();
    for(var i = 0; i < k; i++){
        point_on_curve.x[i] <== tmp[0][i];
        point_on_curve.y[i] <== tmp[1][i];
    }

    component x_check_in_range = CheckInRangeSecp256k1();
    component y_check_in_range = CheckInRangeSecp256k1();
    for(var i = 0; i < k; i++){
        x_check_in_range.in[i] <== tmp[0][i];
        y_check_in_range.in[i] <== tmp[1][i];
    }

    component x3_eq_x1 = BigIsEqual(4);
    for(var i = 0; i < k; i++){
        x3_eq_x1.in[0][i] <== tmp[0][i];
        x3_eq_x1.in[1][i] <== x1[i];
    }
    x3_eq_x1.out === 0;

	for(var i = 0; i < k;i++){
        out[0][i] <== tmp[0][i]*NotZero;
        out[1][i] <== tmp[1][i]*NotZero;
    }

}

template ChosenAddPoint(){
    signal input in[4];
    signal input x[16][4];
    signal input y[16][4];
    signal output px[4];
    signal output py[4];
    signal output IsInfinity;
    
    
    signal v1[2];
    v1[0] <== 1 - in[0];
    v1[1] <== in[0];
    signal v2[4];
    signal v3[8];
    signal v4[16];
    
    for(var j=0;j< 2;j++ ){
        v2[ 2  + j ] <== v1[j] * in[1];
        v2[j] <== v1[j] * (1- in[1]);
    }
    for(var j=0;j< 4;j++ ){
        v3[ 4  + j ] <== v2[j] * in[2];
        v3[j] <== v2[j] * (1- in[2]);
    }
    for(var j=0;j< 8;j++ ){
        v4[ 8  + j ] <== v3[j] * in[3];
        v4[j] <== v3[j] * (1- in[3]);
    }

    signal pxt[16][4];
    for(var i=0;i<16;i++){
        for(var j=0;j<4;j++)
            pxt[i][j] <== v4[i] * x[i][j];
        
    }
    for(var j=0;j<4;j++){
        px[j] <== pxt[0][j] + pxt[1][j] + pxt[2][j] + pxt[3][j] 
                + pxt[4][j] + pxt[5][j] + pxt[6][j] + pxt[7][j]
                + pxt[8][j] + pxt[9][j] + pxt[10][j] + pxt[11][j] 
                + pxt[12][j] + pxt[13][j] + pxt[14][j] + pxt[15][j];
    }

    signal pyt[16][4];
    for(var i=0;i<16;i++){
        for(var j=0;j<4;j++)
            pyt[i][j] <== v4[i] * y[i][j];
    }
    for(var j=0;j<4;j++){
        py[j] <== pyt[0][j] + pyt[1][j] + pyt[2][j] + pyt[3][j] 
                + pyt[4][j] + pyt[5][j] + pyt[6][j] + pyt[7][j]
                + pyt[8][j] + pyt[9][j] + pyt[10][j] + pyt[11][j] 
                + pyt[12][j] + pyt[13][j] + pyt[14][j] + pyt[15][j];
    }

    IsInfinity <== v4[0];  //should be 1 when preadd point is infinity, 0 otherwise
    
}



//use NAF to reduce ec add.
//window = 4;
template Secp256k1ScalarNewMult(n, k) {
    assert(n == 64 && k == 4);
    signal input scalar[k];
    signal input point[2][k];

    signal output out[2][k];

    component n2b[k];
    for (var i = 0; i < k; i++) {
        n2b[i] = Num2Bits(n);
        n2b[i].in <== scalar[i];
    }
    
    signal wx[16][k];
    signal wy[16][k];
    wx[0] <== point[0];   //wx[0] should be equal to point in this part but in further computation wx[0]=0 is used
    wy[0] <== point[1];   
    wx[1] <== point[0];   
    wy[1] <== point[1];
    component preadd[16];
    for(var i=2;i<16;i++){
        preadd[i] = Secp256k1AddUnequal(n,k);
        preadd[i].a[0] <== wx[i-1]; 
        preadd[i].a[1] <== wy[i-1];
        preadd[i].b <== point;
        preadd[i].out[0] ==> wx[i];
        preadd[i].out[1] ==> wy[i]; 
    }


    component chosenPointComp[16][4];
    signal chosenPointX[16][4][4];
    signal chosenPointY[16][4][4];
    signal chosenPointInf[16][4];
    for(var i=15;i>-1;i--){
        for(var j=3;j>-1;j--){
            chosenPointComp[i][j] = ChosenAddPoint();
            chosenPointComp[i][j].in[3] <== n2b[j].out[4*i + 3];
            chosenPointComp[i][j].in[2] <== n2b[j].out[4*i + 2];
            chosenPointComp[i][j].in[1] <== n2b[j].out[4*i + 1];
            chosenPointComp[i][j].in[0] <== n2b[j].out[4*i + 0];
            chosenPointComp[i][j].x <== wx;
            chosenPointComp[i][j].y <== wy;
            chosenPointComp[i][j].px ==> chosenPointX[i][j];
            chosenPointComp[i][j].py ==> chosenPointY[i][j];
            chosenPointComp[i][j].IsInfinity ==> chosenPointInf[i][j];
        }
    }

    component ecdouble[64][4];
    component ecadd[16][4];
    // 257 inter point + 64 + 64
    signal doublepointx[321][4];
    signal doublepointy[321][4];
    signal intermed[2][321][4];
    
    doublepointx[320] <== [0, 0, 0, 0];   //modified double template is used which uses points as input
    doublepointy[320] <== [0, 0, 0, 0];   //

    for(var j=3;j>0;j--){
        for(var i=15;i>-1;i--){
            ecadd[i][j] = Secp256k1AddUnequal(n,k);
            ecadd[i][j].a[0] <== doublepointx[j*80 + i*5 + 5];
            ecadd[i][j].a[1] <== doublepointy[j*80 + i*5 + 5];
            ecadd[i][j].b[0] <== chosenPointX[i][j];
            ecadd[i][j].b[1] <== chosenPointY[i][j];
            ecadd[i][j].out[0] ==> intermed[0][j*80 + i*5 +4];
            ecadd[i][j].out[1] ==> intermed[1][j*80 + i*5 +4];

            for (var j1=0; j1<4; j1++)  // Standard addition
            {
            doublepointx[j*80 + i*5 +4][j1] <== (intermed[0][j*80 + i*5 +4][j1]-doublepointx[j*80 + i*5 + 5][j1])*(1-chosenPointInf[i][j])+doublepointx[j*80 + i*5 + 5][j1];
            doublepointy[j*80 + i*5 +4][j1] <== (intermed[1][j*80 + i*5 +4][j1]-doublepointy[j*80 + i*5 + 5][j1])*(1-chosenPointInf[i][j])+doublepointy[j*80 + i*5 + 5][j1];
            }

            
            //Double is based on Secp256k1Double1(n,k) 
            ecdouble[4*i+3][j] = Secp256k1Double1(n,k);
            ecdouble[4*i+3][j].in[0] <== doublepointx[j*80 + i*5 +4];
            ecdouble[4*i+3][j].in[1] <== doublepointy[j*80 + i*5 +4];
	    ecdouble[4*i+3][j].inbackup[0] <== point[0];
            ecdouble[4*i+3][j].inbackup[1] <== point[1];
            ecdouble[4*i+3][j].out[0] ==> doublepointx[j*80 + i*5 +3];
            ecdouble[4*i+3][j].out[1] ==> doublepointy[j*80 + i*5 +3];


            ecdouble[4*i+2][j] = Secp256k1Double1(n,k);
            ecdouble[4*i+2][j].in[0] <== doublepointx[j*80 + i*5 +3];
            ecdouble[4*i+2][j].in[1] <== doublepointy[j*80 + i*5 +3];
	    ecdouble[4*i+2][j].inbackup[0] <== point[0];
            ecdouble[4*i+2][j].inbackup[1] <== point[1];
            ecdouble[4*i+2][j].out[0] ==> doublepointx[j*80 + i*5 +2];
            ecdouble[4*i+2][j].out[1] ==> doublepointy[j*80 + i*5 +2];


            ecdouble[4*i+1][j] = Secp256k1Double1(n,k);
            ecdouble[4*i+1][j].in[0] <== doublepointx[j*80 + i*5 +2];
            ecdouble[4*i+1][j].in[1] <== doublepointy[j*80 + i*5 +2];
	    ecdouble[4*i+1][j].inbackup[0] <== point[0];
            ecdouble[4*i+1][j].inbackup[1] <== point[1];
            ecdouble[4*i+1][j].out[0] ==> doublepointx[j*80 + i*5 +1];
            ecdouble[4*i+1][j].out[1] ==> doublepointy[j*80 + i*5 +1];


            ecdouble[4*i+0][j] = Secp256k1Double1(n,k);
            ecdouble[4*i+0][j].in[0] <== doublepointx[j*80 + i*5 +1];
            ecdouble[4*i+0][j].in[1] <== doublepointy[j*80 + i*5 +1];
	    ecdouble[4*i+0][j].inbackup[0] <== point[0];
            ecdouble[4*i+0][j].inbackup[1] <== point[1];
            ecdouble[4*i+0][j].out[0] ==> doublepointx[j*80 + i*5 +0];
            ecdouble[4*i+0][j].out[1] ==> doublepointy[j*80 + i*5 +0];

        }
        
    }


    for(var i=15;i>0;i--){
        ecadd[i][0] = Secp256k1AddUnequal(n,k);
        ecadd[i][0].a[0] <== doublepointx[0*80 + i*5 + 5];
        ecadd[i][0].a[1] <== doublepointy[0*80 + i*5 + 5];
        ecadd[i][0].b[0] <== chosenPointX[i][0];
        ecadd[i][0].b[1] <== chosenPointY[i][0];
        ecadd[i][0].out[0] ==> intermed[0][0*80 + i*5 +4];
        ecadd[i][0].out[1] ==> intermed[1][0*80 + i*5 +4];

        for (var j1=0; j1<4; j1++)  //standard additin
            {
            doublepointx[0*80 + i*5 +4][j1] <== (intermed[0][0*80 + i*5 +4][j1]-doublepointx[0*80 + i*5 + 5][j1])*(1-chosenPointInf[i][0])+doublepointx[0*80 + i*5 + 5][j1];
            doublepointy[0*80 + i*5 +4][j1] <== (intermed[1][0*80 + i*5 +4][j1]-doublepointy[0*80 + i*5 + 5][j1])*(1-chosenPointInf[i][0])+doublepointy[0*80 + i*5 + 5][j1];
            }
        
        
       ecdouble[4*i+3][0] = Secp256k1Double1(n,k);
       ecdouble[4*i+3][0].in[0] <== doublepointx[0*80 + i*5 +4];
       ecdouble[4*i+3][0].in[1] <== doublepointy[0*80 + i*5 +4];
       ecdouble[4*i+3][0].inbackup[0] <== point[0];
       ecdouble[4*i+3][0].inbackup[1] <== point[1];
       ecdouble[4*i+3][0].out[0] ==> doublepointx[0*80 + i*5 +3];
       ecdouble[4*i+3][0].out[1] ==> doublepointy[0*80 + i*5 +3];
       

        ecdouble[4*i+2][0] = Secp256k1Double1(n,k);
        ecdouble[4*i+2][0].in[0] <== doublepointx[0*80 + i*5 +3];
        ecdouble[4*i+2][0].in[1] <== doublepointy[0*80 + i*5 +3];
	ecdouble[4*i+2][0].inbackup[0] <== point[0];
        ecdouble[4*i+2][0].inbackup[1] <== point[1];
        ecdouble[4*i+2][0].out[0] ==> doublepointx[0*80 + i*5 +2];
        ecdouble[4*i+2][0].out[1] ==> doublepointy[0*80 + i*5 +2];
        

        ecdouble[4*i+1][0] = Secp256k1Double1(n,k);
        ecdouble[4*i+1][0].in[0] <== doublepointx[0*80 + i*5 +2];
        ecdouble[4*i+1][0].in[1] <== doublepointy[0*80 + i*5 +2];
	ecdouble[4*i+1][0].inbackup[0] <== point[0];
        ecdouble[4*i+1][0].inbackup[1] <== point[1];
        ecdouble[4*i+1][0].out[0] ==> doublepointx[0*80 + i*5 +1];
        ecdouble[4*i+1][0].out[1] ==> doublepointy[0*80 + i*5 +1];


        ecdouble[4*i+0][0] = Secp256k1Double1(n,k);
        ecdouble[4*i+0][0].in[0] <== doublepointx[0*80 + i*5 +1];
        ecdouble[4*i+0][0].in[1] <== doublepointy[0*80 + i*5 +1];
	ecdouble[4*i+0][0].inbackup[0] <== point[0];
        ecdouble[4*i+0][0].inbackup[1] <== point[1];
        ecdouble[4*i+0][0].out[0] ==> doublepointx[0*80 + i*5 +0];
        ecdouble[4*i+0][0].out[1] ==> doublepointy[0*80 + i*5 +0];

    }
    
    ecadd[0][0] = Secp256k1AddUnequal(n,k);
    ecadd[0][0].a[0] <== doublepointx[5];
    ecadd[0][0].a[1] <== doublepointy[5];
    ecadd[0][0].b[0] <== chosenPointX[0][0];
    ecadd[0][0].b[1] <== chosenPointY[0][0];
    ecadd[0][0].out[0] ==> intermed[0][0];
    ecadd[0][0].out[1] ==> intermed[1][0];

    for (var j1=0; j1<4; j1++)  //ovde je obicno sabiranje
            {
            out[0][j1] <== (intermed[0][0][j1]-doublepointx[5][j1])*(1-chosenPointInf[0][0])+doublepointx[5][j1];
            out[1][j1] <== (intermed[1][0][j1]-doublepointy[5][j1])*(1-chosenPointInf[0][0])+doublepointy[5][j1];
            }

}
