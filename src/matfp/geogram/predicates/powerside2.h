/* Automatically generated code, do not edit */
/* Generated from source file: powerside2.pck */

inline int side2_3d_filter( const double* p0, double w0, const double* p1, double w1, const double* p2, double w2, const double* q0, const double* q1) {
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
    double w0_w1 = (w0 - w1);
    double R1;
    R1 = (1 * w0_w1);
    double w0_w2 = (w0 - w2);
    double R2;
    R2 = (1 * w0_w2);
    double L1;
    L1 = (l1 + R1);
    double L2;
    L2 = (l2 + R2);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double a10;
    a10 = (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double a11;
    a11 = (2 * (((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)));
    double a20;
    a20 = (2 * (((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)));
    double a21;
    a21 = (2 * (((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - L1);
    double DeltaLambda1;
    DeltaLambda1 = (L1 - a10);
    double r;
    r = (((Delta * L2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(p1_0_p0_0);
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    } 
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    } 
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    } 
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    } 
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    } 
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    } 
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    } 
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    } 
    else 
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        } 
    } 
    if( (lower_bound_1 < 2.23755023300058943229e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    } 
    else 
    {
        if( (upper_bound_1 > 3.74144419156711063983e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        } 
        eps = (4.44425370757048798480e-15 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        } 
        else 
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            } 
            else 
            {
                return FPG_UNCERTAIN_VALUE;
            } 
        } 
    } 
    Delta_sign = int_tmp_result;
    double max3 = max2;
    if( (max3 < fabs(p2_1_p0_1)) )
    {
        max3 = fabs(p2_1_p0_1);
    } 
    if( (max3 < fabs(p2_2_p0_2)) )
    {
        max3 = fabs(p2_2_p0_2);
    } 
    if( (max3 < fabs(p2_0_p0_0)) )
    {
        max3 = fabs(p2_0_p0_0);
    } 
    double max4 = fabs(w0_w1);
    if( (max4 < fabs(w0_w2)) )
    {
        max4 = fabs(w0_w2);
    } 
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max3;
    upper_bound_1 = max3;
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
    if( ((lower_bound_1 < 1.15014814459393478385e-59) || (max4 < 1.32284075451287083304e-118)) )
    {
        return FPG_UNCERTAIN_VALUE;
    } 
    else 
    {
        if( ((upper_bound_1 > 3.74144419156711063983e+50) || (max4 > 1.39984046386112700994e+101)) )
        {
            return FPG_UNCERTAIN_VALUE;
        } 
        eps = (1.10554268557689364446e-13 * (((max1 * max3) * max3) * std::max( max3, std::max( max1, max4 ) )));
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
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
} 

