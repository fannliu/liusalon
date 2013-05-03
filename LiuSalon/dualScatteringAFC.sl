class dual_scattering_AFC(
	uniform color PrimaryHL_Color        = color(0.86, 0.67, 0.21);
	uniform float PrimaryHL_Intensity    = 0.1;
	uniform float PrimaryHL_LongituShift = -4.5;   //[-10, -5]
	uniform float PrimaryHL_LongituWidth = 0.7;    //[  5, 10]

	uniform color BacklitRim_Color          = color(0.89, 0.98, 0.35);
	uniform float BacklitRim_Intensity      = 0.05; 
	uniform float BacklitRim_LongituShift   = 1;   //-PrimaryHL_LongituShift/2
	uniform float BacklitRim_LongituWidth   = 2;   //PrimaryHL_LongituWidth/2
	uniform float BacklitRim_AzimuthalWidth = 30;

	uniform color SecondaryHL_Color        = color(0.78, 0.4, 0.86);
	uniform float SecondaryHL_Intensity    = 0.1;
	uniform float SecondaryHL_LongituShift = -2;   //-3*PrimaryHL_LongituShift/2
	uniform float SecondaryHL_LongituWidth = 1.4;  //2*PrimaryHL_LongituWidth
	
	uniform float Glints_Intensity = 0.2;          //limit 0.5
	uniform float Glints_AzimuthalShift = 35;      //random per strand[30, 40]
	uniform float Glints_AzimuthalWidth = 0.3;     //[0,1] eqv to frequency

	uniform color ForwardScattering_Color      = color(1,0,0);
	uniform float ForwardScattering_Intensity  = 0.1;
	
	uniform color BackScattering_Color         = color(0,1,0);
	uniform float BackScattering_Intensity     = 0.1;
	uniform float BackScattering_LongituShift  = -2;
	uniform float BackScattering_LongituWidth  = 0.5;
	)
{
	//Defining the uniform parameters used in the pre-computation step
	float R   = 0;
	float TT  = 1;
	float TRT = 2;
	float d_b = 0.7;//backward scattering density factor
	float d_f = 0.7;//forward scattering density factor

	color a_f, a_b, alpha_f, alpha_b, beta_f, beta_b;
	color A_b, delta_b, sigma_b;

	//Defining the varying parameters which have different values for each shading point
	varying float phi;
	varying float theta;
	varying float theta_h;
	varying float hairs_in_front;
	varying color sigma_f;
	varying color T_f;
	
    float g(float deviation, x;)
    {	//unit-integral zero-mean Gaussian distribution
       return exp( - x*x /( 2*deviation*deviation ) ) / ( deviation * sqrt(2*PI) );
    }

	float N(float component;){
		if( component == R)
			return cos(phi * 0.5);
		
		if( component == TT){
			float gamma_TT = radians(BacklitRim_AzimuthalWidth);
			return  g(gamma_TT, PI - phi);
		}
		if( component == TRT){
		
			float G_angle       = radians(Glints_AzimuthalShift);
			float gamma_G       = radians(Glints_AzimuthalWidth);
			        
			float N_TRT_minus_G = cos(phi * 0.5);
			float N_G           = Glints_Intensity * g(gamma_G, G_angle - phi);
			return N_TRT_minus_G + N_G;
		}
	}

	float NG(float component;){
		if( component == R)
			return cos(phi * 0.5);//need integral
		
		if( component == TT){
			float gamma_TT = radians(BacklitRim_AzimuthalWidth);
			return  g(gamma_TT, PI - phi);//need integral
		}
		if( component == TRT){
		
			float G_angle       = radians(Glints_AzimuthalShift);
			float gamma_G       = radians(Glints_AzimuthalWidth);
			        
			float N_TRT_minus_G = cos(phi * 0.5);
			float N_G           = Glints_Intensity * g(gamma_G, G_angle - phi);
			return N_TRT_minus_G + N_G;//need integral
		}
	}

	color M(float component;)
    {
		float alpha, beta;
		if( component == R){
			alpha = radians(PrimaryHL_LongituShift);
			beta  = radians(PrimaryHL_LongituWidth);
		}else if( component == TT){
			alpha = radians(BacklitRim_LongituShift);
			beta  = radians(BacklitRim_LongituWidth);
		}else if( component == TRT){
			alpha = radians(SecondaryHL_LongituShift);
			beta  = radians(SecondaryHL_LongituWidth);
		}
		return g(beta, theta_h - alpha);
    }

	color MG(float component;)
    {
		float alpha, beta;
		if( component == R){
			alpha = radians(PrimaryHL_LongituShift);
			beta  = radians(PrimaryHL_LongituWidth);
		}else if( component == TT){
			alpha = radians(BacklitRim_LongituShift);
			beta  = radians(BacklitRim_LongituWidth);
		}else if( component == TRT){
			alpha = radians(SecondaryHL_LongituShift);
			beta  = radians(SecondaryHL_LongituWidth);
		}
		return g(beta + sigma_f, theta_h - alpha);
    }

	vector GlobalToLocal(vector gv, x, y, z;)
    {
		//transform global vector gv by matrix [LocalUnitX, LocalUnitY, LocalUnitZ]

		float x_ = gv[0] * x[0] + gv[1] * y[0] + gv[2] * z[0];
		float y_ = gv[0] * x[1] + gv[1] * y[1] + gv[2] * z[1];
		float z_ = gv[0] * x[2] + gv[1] * y[2] + gv[2] * z[2];

		return vector(x_, y_, z_);
    }

	//pre-computing and tabulating the uniform variables in the Constructor
	public void construct()
	{
		
		//average forward/backward attenuation
		a_f;
		a_b;

		//average forward/backward scattering shift
		alpha_f;
		alpha_b;

		//average forward/backward scattering deviation
		beta_f;
		beta_b;

		//average backscattering attenuation
		A_b = a_b * pow(a_f,2) / (1 - a_f*a_f) + pow(a_b,3) * pow(a_f,2)/pow(1-a_f*a_f, 2);
		
		//average longitudinal shift
		delta_b = alpha_b * (1 - 2*a_b*a_b/pow(1-a_f*a_f,2)) 
				+ alpha_f * ( 2*pow(1-a_f*a_f,2)+ 4*pow(a_f,2)*pow(a_b,2) )/pow(1-a_f*a_f, 3);

		//average backscattering deviation
		sigma_b = (1 + d_b * pow(a_f,2)) 
				*( a_b*sqrt(2*pow(beta_f,2) + pow(beta_b,2)) + pow(a_b,3)*sqrt(2*pow(beta_f,2) + pow(beta_b,2)) )
				/ (a_b + pow(a_b,3)*(2*beta_f + 3*beta_b);
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
		float alpha_back = radians(BackScattering_LongituShift);
		float beta_back  = radians(BackScattering_LongituWidth);
		
		illuminance(P) //P is the shading point position, a function of (u,v)
		{
			vector omega_i = GlobalToLocal( normalize(L), lx, ly, lz ); //L is light ray (from shading point to the light source)
			float phi_i    = atan(omega_i[1], omega_i[0]);
				  phi      = abs(phi_o - phi_i); //relative azimuth (within the normal plane)

			if ( phi > PI )
				phi -= 2 * PI;
			phi = abs(phi);

			float theta_i = PI * 0.5 - acos(omega_i[2]);
				  theta   = theta_i + theta_o;
				  theta_h = theta * 0.5; //half longitudial angle (wrt the normal plane)

			//compute the amount of shadow from the deep shadow maps?????
			float shadowed = 0;
			lightsource("out_shadow", shadowed);
			float illuminated = 1 - shadowed;
			//estimate the number of hairs in front of the shading point
			hairs_in_front = shadowed * hairs_that_cast_full_shadow;
			//use the number of hairs in front of the shading point to approximate sigma_f
			sigma_f = hairs_in_front * beta_f;
			//use the number of hairs in front of the shading point to approximate T_f
			T_f = d_f * pow(a_f, hairs_in_front);


			
			//backscattering for direct and indirect lighting
			color f_direct_back  =  2 * A_b * g( sigma_b + beta_back, theta_h - delta_b + alpha_back) 
									/ (PI * cos(theta) * cos(theta));
				  f_direct_back = BackScattering_Color * BackScattering_Intensity * f_direct_back;

			color f_scatter_back = 2 * A_b * g( sigma_b + sigma_f + beta_back, theta_h - delta_b + alpha_back) 
									/ (PI * cos(theta) * cos(theta));
		          f_scatter_back = BackScattering_Color * BackScattering_Intensity * f_scatter_back;

			//single scattering for direct and indirect lighting
			color f_direct_s  = M(R)*N(R) + M(TT)*N(TT) + M(TRT)*N(TRT);
			color f_scatter_s = MG(R)*NG(R) + MG(TT)*NG(TT) + MG(TRT)*NG(TRT);
			      f_scatter_s = ForwardScattering_Color * ForwardScattering_Intensity * f_scatter_s;
			
			
			color F_direct = illuminated * ( f_direct_s + f_direct_back);
			color F_scatter = (T_f - illuminated) * d_f * ( f_scatter_s + PI * d_b * f_scatter_back);

			//combine the direct and indirect scattering components
			Ci  += (F_direct + F_scatter) * cos(theta_i);
		}

		Oi  = Os; //Os is the surface opacity
		Ci *= Oi;
	
    }
}