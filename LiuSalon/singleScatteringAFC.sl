class single_scattering_AFC(
	uniform color colorR = color(0.86, 0.67, 0.21);
	uniform float intensityR = 0.3;
	uniform float longitudinalShiftR = -7.5;//[-10,-5]
	uniform float longitudinalWidthR = 5.0;//[5, 10]

	uniform color colorTT = color(0, 1, 0);
	uniform float intensityTT = 0.5;
	uniform float longitudinalShiftTT = 3.75;
	uniform float longitudinalWidthTT = 3.75;
	uniform float azimuthalWidthTT = 3;

	uniform color colorTRT = color(0, 0, 1);
	uniform float intensityTRT = 0.5;
	uniform float longitudinalShiftTRT = 11.25;
	uniform float longitudinalWidthTRT = 15;
	
	uniform float intensityG = 0.2;
	uniform float azimuthalWidthG = 10;)//equivalent to frequency
{
	//unit-integral zero-mean Gaussian function
    float g(float deviation, x;)
    {
       return exp(-x*x/(2*deviation*deviation))/(deviation*sqrt(2*PI));
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
	vector GlobalToLocal(vector l, x, y, z;)
    {
        float a = l[0] * x[0] + l[1] * y[0] + l[2] * z[0];
        float b = l[0] * x[1] + l[1] * y[1] + l[2] * z[1];
        float c = l[0] * x[2] + l[1] * y[2] + l[2] * z[2];

        return vector(a,b,c);
    }
    public void surface(output color Ci, Oi;)
    {
		// Get unit vectors along local coordinates
		vector lx  =   normalize(dPdu);
		vector ly  =   normalize(N);//the shading normal
		vector lz  =   normalize(dPdv);	
		
        
        vector omega_r = GlobalToLocal(-normalize(I), lx, ly, lz);//I is the incident ray dir(from eye to the shading point)
		
		color singleScatteringResult = 0;

		illuminance(P)//P is the shading point position, a function ofthe surface parameters (u,v)
		{
			vector omega_i = GlobalToLocal(normalize(L), lx, ly, lz);//light ray (from shading point to the light source)
			 
			float phi_i = atan(omega_i[1], omega_i[0]);
			float phi_r = atan(omega_r[1], omega_r[0]);

			float phi = abs(phi_r -  phi_i);//relative azimuth
            if ( phi > PI )
                phi -= 2 * PI;
            phi = abs(phi);

			float theta_i = PI * 0.5 - acos(omega_i[2]);
			float theta_r = PI * 0.5 - acos(omega_r[2]);
			float theta = theta_i + theta_r;
			
			float theta_h =  theta * 0.5;//half angle btw theta_i and theta_o; azimuthal angle (dir wrt the normal plane)

			color f_R = R(theta_h, phi);
			color f_TT = TT(theta_h, phi);
			color f_TRT = TRT(theta_h, phi);
			//singleScatteringResult += (f_R  + f_TT + f_TRT)/(cos(theta)*cos(theta));
			singleScatteringResult += (f_R)/(cos(theta)*cos(theta));
			//singleScatteringResult += (f_TT)/(cos(theta)*cos(theta));
			//singleScatteringResult += (f_TRT)/(cos(theta)*cos(theta));
		
		}
		Oi = Os;//surface opacity
		Ci = singleScatteringResult * Oi;
	
    }
}