class single_scattering_AFC(
	uniform color PrimaryHL_Color        = color(0.95,0.85,0.64);
	uniform float PrimaryHL_Intensity    = 0.4;
	uniform float PrimaryHL_LongituShift = -4.5;   //[-10, -5]
	uniform float PrimaryHL_LongituWidth = 2;    //[  5, 10]

	uniform color BacklitRim_Color          = color(0.80,0.67,0.39);
	uniform float BacklitRim_Intensity      = 0.15; 
	uniform float BacklitRim_LongituShift   = 2.25;   //-PrimaryHL_LongituShift/2
	uniform float BacklitRim_LongituWidth   = 1;   //PrimaryHL_LongituWidth/2
	uniform float BacklitRim_AzimuthalWidth = 30;

	uniform color SecondaryHL_Color        = color(0.78,0.61,0.21);
	uniform float SecondaryHL_Intensity    = 0.3;
	uniform float SecondaryHL_LongituShift = -4;   //-3*PrimaryHL_LongituShift/2
	uniform float SecondaryHL_LongituWidth = 4;  //2*PrimaryHL_LongituWidth
	
	uniform float Glints_Intensity = 0.3;          //limit 0.5
	uniform float Glints_AzimuthalShift = 35;      //random per strand[30, 40]
	uniform float Glints_AzimuthalWidth = 0.3;     //[0,1] eqv to frequency
	uniform float attenuationFromRoot = 1.0;
	 
	uniform float Diffuse_Intensity = 0.2;
	uniform color Diffuse_Color = color(1, 1, 1);
	)
{
	//unit-integral zero-mean Gaussian distribution
    float g(float deviation, x;)
    {
       return exp( - x*x /( 2*deviation*deviation ) );// / ( deviation * sqrt(2*PI) );
    }

    color R(float theta_h, phi;)
    {
		float alpha_R = radians(PrimaryHL_LongituShift);
		float beta_R  = radians(PrimaryHL_LongituWidth);

		float M_R     = g(beta_R, theta_h - alpha_R);
		float N_R     = cos(phi * 0.5);

		return PrimaryHL_Color * PrimaryHL_Intensity * M_R * N_R;
    }

    color TT(float theta_h, phi;)
	{
		float alpha_TT = radians(BacklitRim_LongituShift);
		float beta_TT  = radians(BacklitRim_LongituWidth);
		float gamma_TT = radians(BacklitRim_AzimuthalWidth);

		float M_TT     = g(beta_TT, theta_h - alpha_TT);
		float N_TT     = g(gamma_TT, PI - phi);
    
		return BacklitRim_Color * BacklitRim_Intensity * M_TT * N_TT;
    }

    color TRT(float theta_h, phi;)
    {
		float alpha_TRT     = radians(SecondaryHL_LongituShift);
		float beta_TRT      = radians(SecondaryHL_LongituWidth);

		float G_angle       = radians(Glints_AzimuthalShift);
		float gamma_G       = radians(Glints_AzimuthalWidth);

		float M_TRT         = g(beta_TRT, theta_h - alpha_TRT);
        
		float N_TRT_minus_G = cos(phi * 0.5);
		float N_G           = Glints_Intensity * g(gamma_G, G_angle - phi);
		float N_TRT         = N_TRT_minus_G + N_G;
        
		return SecondaryHL_Intensity * SecondaryHL_Color * M_TRT * N_TRT;
    }

	vector GlobalToLocal(vector gv, x, y, z;)
    {
		//transform global vector gv by matrix [LocalUnitX, LocalUnitY, LocalUnitZ]

		float x_ = gv[0] * x[0] + gv[1] * y[0] + gv[2] * z[0];
		float y_ = gv[0] * x[1] + gv[1] * y[1] + gv[2] * z[1];
		float z_ = gv[0] * x[2] + gv[1] * y[2] + gv[2] * z[2];

		return vector(x_, y_, z_);
    }

    public void surface(output color Ci, Oi;)
    {
	
		// Get unit vectors along local axis in globle coordinate system
		vector lx  =  normalize(dPdu);
		vector ly  =  normalize(   N);  //the shading normal
		vector lz  =  normalize(dPdv);	//hair tangent (from root to tip)
		
        
		vector omega_o = GlobalToLocal( -normalize(I), lx, ly, lz ); //I is the incident ray dir(from eye to the shading point)
		float    phi_o = atan(omega_o[1], omega_o[0]);
		float  theta_o = PI * 0.5 - acos(omega_o[2]);

		color singleScatteringResult = 0;

		color diffuseResult = 0;
		illuminance(P) //P is the shading point position, a function of (u,v)
		{
			vector Ln = normalize(L);
			vector omega_i = GlobalToLocal( Ln, lx, ly, lz ); //L is light ray (from shading point to the light source)
			float phi_i    = atan(omega_i[1], omega_i[0]);
			float phi      = abs(phi_o - phi_i); //relative azimuth (within the normal plane)

			if ( phi > PI )
				phi -= 2 * PI;
			phi = abs(phi);

			float theta_i = PI * 0.5 - acos(omega_i[2]);
			float theta   = theta_i + theta_o;
			float theta_h = theta * 0.5; //half longitudial angle (wrt the normal plane)
			
			
			float cosTL = clamp(Ln.ly, 0.0, 1.0);
			float sinTL = sqrt(1 - cosTL * cosTL);
			diffuseResult += Cl * Diffuse_Intensity * Diffuse_Color * sinTL;
			
			color f_R   =   R(theta_h, phi);
			color f_TT  =  TT(theta_h, phi);
			color f_TRT = TRT(theta_h, phi);
			singleScatteringResult += ( f_R  + f_TT + f_TRT ) / ( cos(theta) * cos(theta) );
			
			 
		}

		Oi = Os; //Os is the surface opacity
		Ci = (diffuseResult * pow(v, attenuationFromRoot) + singleScatteringResult) * Oi;
	
    }
}