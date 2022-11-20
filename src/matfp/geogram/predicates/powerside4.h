/* Automatically generated code, do not edit */
/* Generated from source file: powerside4.pck */

inline int side4_3d_filter( const double* p0, double w0, const double* p1, double w1, const double* p2, double w2, const double* p3, double w3, const double* p4, double w4) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double l1;
    l1 = (1 * (((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double l2;
    l2 = (1 * (((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double l3;
    l3 = (1 * (((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double l4;
    l4 = (1 * (((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)));
    double w0_w1 = (w0 - w1);
    double R1;
    R1 = (1 * w0_w1);
    double w0_w2 = (w0 - w2);
    double R2;
    R2 = (1 * w0_w2);
    double w0_w3 = (w0 - w3);
    double R3;
    R3 = (1 * w0_w3);
    double w0_w4 = (w0 - w4);
    double R4;
    R4 = (1 * w0_w4);
    double L1;
    L1 = (l1 + R1);
    double L2;
    L2 = (l2 + R2);
    double L3;
    L3 = (l3 + R3);
    double L4;
    L4 = (l4 + R4);
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a14;
    a14 = -L1;
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double a24;
    a24 = -L2;
    double a31;
    a31 = (p3[0] - p0[0]);
    double a32;
    a32 = (p3[1] - p0[1]);
    double a33;
    a33 = (p3[2] - p0[2]);
    double a34;
    a34 = -L3;
    double a41;
    a41 = (p4[0] - p0[0]);
    double a42;
    a42 = (p4[1] - p0[1]);
    double a43;
    a43 = (p4[2] - p0[2]);
    double a44;
    a44 = -L4;
    double Delta1;
    Delta1 = (((a21 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a22 * a43) - (a23 * a42)))) + (a41 * ((a22 * a33) - (a23 * a32))));
    double Delta2;
    Delta2 = (((a11 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a12 * a43) - (a13 * a42)))) + (a41 * ((a12 * a33) - (a13 * a32))));
    double Delta3;
    Delta3 = (((a11 * ((a22 * a43) - (a23 * a42))) - (a21 * ((a12 * a43) - (a13 * a42)))) + (a41 * ((a12 * a23) - (a13 * a22))));
    double Delta4;
    Delta4 = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    double r;
    r = ((((Delta1 * a14) - (Delta2 * a24)) + (Delta3 * a34)) - (Delta4 * a44));
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a21)) )
    {
        max1 = fabs(a21);
    } 
    if( (max1 < fabs(a31)) )
    {
        max1 = fabs(a31);
    } 
    double max2 = fabs(a12);
    if( (max2 < fabs(a13)) )
    {
        max2 = fabs(a13);
    } 
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    } 
    if( (max2 < fabs(a23)) )
    {
        max2 = fabs(a23);
    } 
    double max3 = fabs(a22);
    if( (max3 < fabs(a23)) )
    {
        max3 = fabs(a23);
    } 
    if( (max3 < fabs(a32)) )
    {
        max3 = fabs(a32);
    } 
    if( (max3 < fabs(a33)) )
    {
        max3 = fabs(a33);
    } 
    double lower_bound_1;
    double upper_bound_1;
    int Delta4_sign;
    int int_tmp_result;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    } 
    else 
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        } 
    } 
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    } 
    else 
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        } 
    } 
    if( (lower_bound_1 < 1.63288018496748314939e-98) )
    {
        return FPG_UNCERTAIN_VALUE;
    } 
    else 
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        } 
        eps = (5.11071278299732992696e-15 * ((max2 * max3) * max1));
        if( (Delta4 > eps) )
        {
            int_tmp_result = 1;
        } 
        else 
        {
            if( (Delta4 < -eps) )
            {
                int_tmp_result = -1;
            } 
            else 
            {
                return FPG_UNCERTAIN_VALUE;
            } 
        } 
    } 
    Delta4_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max4 = max1;
    if( (max4 < fabs(a41)) )
    {
        max4 = fabs(a41);
    } 
    double max5 = max2;
    if( (max5 < max3) )
    {
        max5 = max3;
    } 
    double max6 = max3;
    if( (max6 < fabs(a42)) )
    {
        max6 = fabs(a42);
    } 
    if( (max6 < fabs(a43)) )
    {
        max6 = fabs(a43);
    } 
    double max7 = fabs(p1_2_p0_2);
    if( (max7 < fabs(p3_2_p0_2)) )
    {
        max7 = fabs(p3_2_p0_2);
    } 
    if( (max7 < fabs(p2_2_p0_2)) )
    {
        max7 = fabs(p2_2_p0_2);
    } 
    if( (max7 < fabs(p4_0_p0_0)) )
    {
        max7 = fabs(p4_0_p0_0);
    } 
    if( (max7 < fabs(p2_1_p0_1)) )
    {
        max7 = fabs(p2_1_p0_1);
    } 
    if( (max7 < fabs(p4_1_p0_1)) )
    {
        max7 = fabs(p4_1_p0_1);
    } 
    if( (max7 < fabs(p3_1_p0_1)) )
    {
        max7 = fabs(p3_1_p0_1);
    } 
    if( (max7 < fabs(p2_0_p0_0)) )
    {
        max7 = fabs(p2_0_p0_0);
    } 
    if( (max7 < fabs(p1_1_p0_1)) )
    {
        max7 = fabs(p1_1_p0_1);
    } 
    if( (max7 < fabs(p3_0_p0_0)) )
    {
        max7 = fabs(p3_0_p0_0);
    } 
    if( (max7 < fabs(p4_2_p0_2)) )
    {
        max7 = fabs(p4_2_p0_2);
    } 
    if( (max7 < fabs(p1_0_p0_0)) )
    {
        max7 = fabs(p1_0_p0_0);
    } 
    double max8 = fabs(w0_w1);
    if( (max8 < fabs(w0_w2)) )
    {
        max8 = fabs(w0_w2);
    } 
    if( (max8 < fabs(w0_w3)) )
    {
        max8 = fabs(w0_w3);
    } 
    if( (max8 < fabs(w0_w4)) )
    {
        max8 = fabs(w0_w4);
    } 
    lower_bound_1 = max7;
    upper_bound_1 = max7;
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    } 
    else 
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        } 
    } 
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    } 
    else 
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        } 
    } 
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    } 
    else 
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        } 
    } 
    if( ((lower_bound_1 < 7.07260400778147882026e-50) || (max8 < 5.00217274508866438928e-99)) )
    {
        return FPG_UNCERTAIN_VALUE;
    } 
    else 
    {
        if( ((upper_bound_1 > 1.87072209578355531992e+50) || (max8 > 3.49960115965281752486e+100)) )
        {
            return FPG_UNCERTAIN_VALUE;
        } 
        eps = (1.77774053336465084468e-13 * ((((max5 * max6) * max4) * max7) * std::max( max7, max8 )));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        } 
        else 
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            } 
            else 
            {
                return FPG_UNCERTAIN_VALUE;
            } 
        } 
    } 
    return (Delta4_sign * int_tmp_result_FFWKCAA);
} 

