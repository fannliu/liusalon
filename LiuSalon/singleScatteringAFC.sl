class single_scattering_AFC(
				 uniform color colorR = color(1, 1, 1);
                 uniform float intensityR = 1;
                 uniform float longitudinalShiftR = -7.5;
                 uniform float longitudinalWidthR = 7.5;

                 uniform color colorTT = color(1, 1, 1);
				 uniform float intensityTT = 1;
                 uniform float longitudinalShiftTT = 3.75;
                 uniform float longitudinalWidthTT = 3.75;
                 uniform float azimuthalWidthTT = 3;

				 uniform color colorTRT = color(1, 1, 1);
                 uniform float intensityTRT = 1;
                 uniform float longitudinalShiftTRT = 11.25;
                 uniform float longitudinalWidthTRT = 15;

				 uniform float intensityG = 1;
				 uniform float azimuthalWidthG = 10;//equivalent to frequency
				 uniform float azimuthalShiftG = 40;//Gangle????different per strand
{
    float UnitGaussian(float sigma, x_minus_mu)
    {
        return exp(- pow(x_minus_mu, 2.0) * 0.5 / pow(sigma, 2.0));
    }

    color R(float theta_h, longitudinalShiftR, longitudinalWidthR, phi)
    {
        float M_R = UnitGaussian(longitudinalWidthR, theta_h - longitudinalShiftR);
        float N_R = cos(phi * 0.5);
        return colorR * intensityR * M_R * N_R;
    }

    color TT(float theta_h, longitudinalShiftTT, longitudinalWidthTT, phi, azimuthalWidthTT)
    {
        float M_TT = UnitGaussian(longitudinalWidthTT, theta_h - longitudinalShiftTT);
        float N_TT = UnitGaussian(azimuthalWidthTT, PI - phi);
    
        return colorTT * intensityTT * M_TT * N_TT;
    }

    color TRT(float theta_h, longitudinalShiftTRT, longitudinalWidthTRT, phi, azimuthalShiftG, azimuthalWidthG)
    {
        float M_TRT = UnitGaussian(longitudinalWidthTRT, theta_h - longitudinalShiftTRT);
        
		float N_TRT_minus_G = cos(phi * 0.5);
        float N_G = intensityG * UnitGaussian(azimuthalWidthG, azimuthalShiftG - phi);
        float N_TRT = N_TRT_minus_G + N_G;
        
        return intensityTRT * colorTRT * M_TRT * N_TRT;
    }

    public void surface(output color Ci, Oi)
    {

    }
}