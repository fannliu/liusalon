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
	vector GlobalToLocal(vector l, x, y, z;)
    {
        float a = l[0] * x[0] + l[1] * y[0] + l[2] * z[0];
        float b = l[0] * x[1] + l[1] * y[1] + l[2] * z[1];
        float c = l[0] * x[2] + l[1] * y[2] + l[2] * z[2];

        return vector(a,b,c);
    }

	float angleBtwVec(vector vi, vo; int axis)
	{
		vi[axis] = 0;
		vo[axis] = 0;

		float angle = acos(vi * vo);

		//clamp angle between [-Pi, Pi]
		return angle;
	}
    public void surface(output color Ci, Oi;)
    {
		// Get local frame
		vector U  =   normalize(dPdu);
        vector Nn =   normalize(N);//the shading normal
		vector V  =	  normalize(dPdv);
        
        vector omega_r = GlobalToLocal(-normalize(I), U, Nn, V);//I is the incident ray dir(from eye to the shading point) in local coordinate
		
		color singleScatteringResult = 0;

		illuminance(P)//P is the shading point position, a function ofthe surface parameters (u,v)
		{
			vector omega_i = GlobalToLocal(normalize(L), U, Nn, V);//light ray (from shading point to the light source)
			

			float phi =  angleBtwVec(omega_i,omega_o, 2);//longitudinal inclination (dir within the normal plane)
			float theta_h =  angleBtwVec(omega_i, omega_o, 0);//half angle btw theta_i and theta_o; azimuthal angle (dir wrt the normal plane)

			color f_R = R(theta_h, phi);
			color f_TT = TT(theta_h, phi);
			color f_TRT = TRT(theta_h, phi);
			singleScatteringResult += (f_R  + f_TT + f_TRT)/(cos(theta)*cos(theta));
		}
		Oi = Os;
		Ci = singleScatteringResult * Oi;
	
    }
};