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
{
    float g(float variance, x;)
    {
        return exp(-x*x*0.5/variance)/sqrt(2*PI*variance);
    }

    color R(float theta_h, phi;)
    {
		float alpha_R = radians(longitudinalShiftR);
        float beta_R  = radians(longitudinalWidthR);

        float M_R = g(beta_R, theta_h - alpha_R);
        float N_R = cos(phi * 0.5);
        return colorR * intensityR * M_R * N_R;
    }

    color TT(float theta_h, phi;)
	{

		float alpha_TT = radians(longitudinalShiftTT);
        float beta_TT  = radians(longitudinalWidthTT);
        float gamma_TT = radians(azimuthalWidthTT);

        float M_TT = g(beta_TT, theta_h - beta_TT);
        float N_TT = g(gamma_TT, PI - phi);
    
        return colorTT * intensityTT * M_TT * N_TT;
    }

    color TRT(float theta_h, phi;)
    {
		float alpha_TRT = radians(longitudinalShiftTRT);
        float beta_TRT  = radians(longitudinalWidthTRT);

        float G_angle = radians(30);//random between 30-45
        float gamma_G = radians(azimuthalWidthG);

        float M_TRT = g(beta_TRT, theta_h - alpha_TRT);
        
		float N_TRT_minus_G = cos(phi * 0.5);
        float N_G = intensityG * g(gamma_G, G_angle - phi);
        float N_TRT = N_TRT_minus_G + N_G;
        
        return intensityTRT * colorTRT * M_TRT * N_TRT;
    }
	vector LocalSpherical(vector l, x, y, z)
    {
        float a = l[0] * x[0] + l[1] * y[0] + l[2] * z[0];
        float b = l[0] * x[1] + l[1] * y[1] + l[2] * z[1];
        float c = l[0] * x[2] + l[1] * y[2] + l[2] * z[2];

        return vector(atan(b, a), PI * 0.5 - acos(c), 0);
    }
    public void surface(output color Ci, Oi;)
    {
		
		vector u = -normalize(dPdv);//tangent to the hair fiber

		// Get local frame.-need more work
      
        vector T  =   normalize(dPdv);
        vector Nn =   normalize(N);
        vector B  =   normalize(Nn ^ T);
        vector Vn = - normalize(I);
		
		color singleScatteringResult = 0;

		illuminance(P)
		{
			vector Ln = normalize(L);
			vector omega_i = LocalSpherical(Ln, B, Nn, T);//incoming light direction
			vector omega_o = LocalSpherical(Vn, B, Nn, T);//outgoing view direction
		
			float phi = abs(omega_i[0] - omega_o[0]);//longitudinal inclination (dir w.r.t. to the normal plane)
			if( phi > PI )
			phi -= 2 * PI;//clamp phi to [-pi, pi]
			phi = abs(phi);

			float theta_h = abs(omega_o[1] - omega-o[1]) * 0,5;//half angle btw theta_i and theta_o; azimuthal angle (dir within the normal plane)

			color f_R = R(theta_h, phi);
			color f_TT = TT(theta_h, phi);
			color f_TRT = TRT(theta_h, phi);
			singleScatteringResult += (f_R  + f_TT + f_TRT)/(cos(theta)*cos(theta));
		}
		Oi = Os;
		Ci = singleScatteringResult * Oi;
	
    }
}