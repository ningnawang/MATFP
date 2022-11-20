/* Automatically generated code, do not edit */
/* Generated from source file: powerside1.pck */

inline int side1_power_3d_filter( const double* p0, double w0, const double* p1, double w1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double r;
    r = (1 * (((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)));
    double w0_w1 = (w0 - w1);
    r = (r + w0_w1);
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    r = (r - (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0_2_p1_2);
    if( (max1 < fabs(p0_0_p1_0)) )
    {
        max1 = fabs(p0_0_p1_0);
    } 
    if( (max1 < fabs(p0_1_p1_1)) )
    {
        max1 = fabs(p0_1_p1_1);
    } 
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    } 
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    } 
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    } 
    double max2 = fabs(p0_2_p1_2);
    if( (max2 < fabs(p0_0_p1_0)) )
    {
        max2 = fabs(p0_0_p1_0);
    } 
    if( (max2 < fabs(p0_1_p1_1)) )
    {
        max2 = fabs(p0_1_p1_1);
    } 
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    } 
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    } 
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    } 
    double max3 = fabs(w0_w1);
    double lower_bound_1;
    double upper_bound_1;
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
    if( ((lower_bound_1 < 2.08650753039307467708e-147) || (max3 < 4.35351367438700797978e-294)) )
    {
        return FPG_UNCERTAIN_VALUE;
    } 
    else 
    {
        if( ((upper_bound_1 > 5.59936185544450928309e+101) || (max3 > 3.13528531882069776730e+203)) )
        {
            return FPG_UNCERTAIN_VALUE;
        } 
        eps = (5.11098396588934852270e-15 * std::max( max3, (max1 * max2) ));
        if( (r > eps) )
        {
            int_tmp_result = 1;
        } 
        else 
        {
            if( (r < -eps) )
            {
                int_tmp_result = -1;
            } 
            else 
            {
                return FPG_UNCERTAIN_VALUE;
            } 
        } 
    } 
    return int_tmp_result;
} 

