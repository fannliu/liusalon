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
	constant float R   = 0;
	constant float TT  = 1;
	constant float TRT = 2;
	
	constant float hemi_f = 0;
	constant float hemi_b = 1;
	
	
	constant float M_PI_2 = PI * 0.5;
	constant float M_2_PI = 2.0 / PI;
	constant float segment = 0.1;
	constant float tableSize;
	
	constant float d_b = 0.7;//backward scattering density factor
	constant float d_f = 0.7;//forward scattering density factor
	
	constant float hairs_that_cast_full_shadow = 15;//not sure??????
	
	//Defining the uniform parameters used in the pre-computation step
	uniform color a_f[];
	uniform color a_b[];
	uniform color alpha_f[];
	uniform color alpha_b[];
	uniform color beta_f[];
	uniform color beta_b[];
	
	uniform color A_b[];
	uniform color delta_b[];
	uniform color sigma_b[];

	//Defining the varying parameters which have different values for each shading point
	varying float hairs_in_front;
	varying color sigma_f;
	varying color T_f;
	
    float g(float deviation, x;)
    {	//unit-integral zero-mean Gaussian distribution
       return exp( - x*x /( 2*deviation*deviation ) ) / ( deviation * sqrt(2*PI) );
    }
	color M(uniform float component; float theta_h;)
    {
		float alpha_, beta_;
		if( component == R){
			alpha_ = radians(PrimaryHL_LongituShift);
			beta_  = radians(PrimaryHL_LongituWidth);
		}else if( component == TT){
			alpha_ = radians(BacklitRim_LongituShift);
			beta_  = radians(BacklitRim_LongituWidth);
		}else if( component == TRT){
			alpha_ = radians(SecondaryHL_LongituShift);
			beta_  = radians(SecondaryHL_LongituWidth);
		}
		return g(beta, theta_h - alpha_);
    }

	color MG(uniform float component; float theta_h;)
    {
		float alpha_, beta_;
		if( component == R){
			alpha_ = radians(PrimaryHL_LongituShift);
			beta_  = radians(PrimaryHL_LongituWidth);
		}else if( component == TT){
			alpha_ = radians(BacklitRim_LongituShift);
			beta_  = radians(BacklitRim_LongituWidth);
		}else if( component == TRT){
			alpha_ = radians(SecondaryHL_LongituShift);
			beta_  = radians(SecondaryHL_LongituWidth);
		}
		return g(beta_ + sigma_f, theta_h - alpha_);
    }

	float N_(uniform float component; float phi;){
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

	float NG(uniform float component;){
		float result = 0;
		if( component == R){
			for(float phi = M_PI_2; phi <= PI; phi += segment)
				result += cos(phi * 0.5);
		}
		else if( component == TT){
			float gamma_TT = radians(BacklitRim_AzimuthalWidth);
			for(float phi = M_PI_2; phi <= PI; phi += segment)
				result += g(gamma_TT, PI - phi);
		}
		else if( component == TRT){
		
			float G_angle       = radians(Glints_AzimuthalShift);
			float gamma_G       = radians(Glints_AzimuthalWidth);
			
			for(float phi = M_PI_2; phi <= PI; phi += segment){        
				float N_TRT_minus_G = cos(phi * 0.5);
				float N_G           = Glints_Intensity * g(gamma_G, G_angle - phi);
				result += (N_TRT_minus_G + N_G);
			}
		}
		return result * M_2_PI;
	}

	

	vector GlobalToLocal(vector gv, x, y, z;)
    {
		//transform global vector gv by matrix [LocalUnitX, LocalUnitY, LocalUnitZ]

		float x_ = gv[0] * x[0] + gv[1] * y[0] + gv[2] * z[0];
		float y_ = gv[0] * x[1] + gv[1] * y[1] + gv[2] * z[1];
		float z_ = gv[0] * x[2] + gv[1] * y[2] + gv[2] * z[2];

		return vector(x_, y_, z_);
    }
	color singleScatteringIntegral()
	{
	}
	color singleScattering(float theta, phi;)
	{
		float theta_h = theta * 0.5;
		color f_R   =   PrimaryHL_Color *   PrimaryHL_Intensity * M(R, theta_h)   * N_(R, phi);
		color f_TT  =  BacklitRim_Color *  BacklitRim_Intensity * M(TT, theta_h)  * N_(TT, phi);
		color f_TRT = SecondaryHL_Color * SecondaryHL_Intensity * M(TRT, theta_h) * N_(TRT, phi);

		return (f_R + f_TT + f_TRT)/(cos(theta) * cos(theta))
	}
	
	color normalizedSingleScattering(float theta, phi;)
	{
		return singleScattering(theta, phi) / singleScatteringIntegral();
	}

	//pre-computing and tabulating the uniform variables in the Constructor
	void populate_a(uniform float hemisphere; output uniform color a[]) 
    {

    }
	
	color integrateOverHemisphere(uniform float theta; uniform float hemisphere){
	}
	color integrateOverHemisphereWeighted(uniform float theta; uniform float coef; uniform float hemisphere){
	
	}
	
	void populate_alphabeta(uniform float hemisphere; output uniform color target[]; float coef;)
    {
        uniform float theta = - M_PI_2;
        uniform float i;
        uniform color denominator = 0.0;
        uniform color numerator = 0.0;
 
        for (i = 0; i < tableSize; i++) {
            denominator = integrateOverHemisphere(theta, hemisphere);
              numerator = integrateOverHemisphereWeighted(theta, coef, hemisphere);
 
            push( target, numerator / denominator );
            theta += segment;
        }
    }

	void populate_A_b(output uniform color A[])
    {
        uniform float i;
 
        for (i = 0; i < tableSize; i += 1) {
            uniform color af = a_f[i];
            uniform color ab = a_b[i];
			
			uniform color afPow2 = pow(af, 2);
			
			Ab = ab * afPow2 / (1 - afPow2) + pow(ab,3) * afPow2/pow(1-afPow2, 2);
			
            push( A, Ab );
        }
    }
	
	void populate_delta_b(output uniform color delta[]) 
    {
        uniform float i;
 
        for (i = 0; i < tableSize; i++) {
            uniform color af = a_f[i];
            uniform color ab = a_b[i];
			
			uniform color alphaf = alpha_f[i];
            uniform color alphab = alpha_b[i];
			
			uniform color afPow2 = pow(af, 2);
			uniform color abPow2 = pow(ab, 2);
			
			uniform color deltab = alphab * (1 - 2*abPow2/pow(1-afPow2,2)) 
				+ alphaf * ( 2*pow(1-afPow2,2)+ 4*afPow2*abPow2 ) / pow(1-afPow2, 3);
 
            push( delta, deltab );
        }
    }
	
	void populate_sigma_b(output uniform color sigma[])
    {
        uniform float i;
 
        for (i = 0; i < tableSize; i++) {
            uniform color af = a_f[i];
            uniform color ab = a_b[i];
            
            uniform color betaf = beta_f[i];
            uniform color betab = beta_b[i];
			
			uniform color betabPow2 = pow(betab, 2);
			uniform color betafPow2 = pow(betaf, 2);
			uniform color abPow3 = pow(ab,3);
			
			uniform color sigmab = (1 + d_b * pow(af,2)) * (ab + abPow3) * sqrt(2*betafPow2 + betabPow2)
									/ ( ab + abPow3*(2*betaf + 3*betab) );
				
            push( sigma, sigmab );
        }
    }
	
	public void construct()
	{	
	
		tableSize = ceil(PI/segment);
		
		/*average forward/backward attenuation*/
		reserve(a_f, tableSize);
		reserve(a_b, tableSize);
		populate_a(hemi_f, a_f);
		populate_a(hemi_b, a_b);

		/*average forward/backward scattering shift*/
		reserve(alpha_f, tableSize);
		reserve(alpha_b, tableSize);
		populate_alphabeta(hemi_f, alpha_f, 0);
		populate_alphabeta(hemi_b, alpha_b, 0);

		/*average forward/backward scattering deviation/width*/
		reserve(beta_f, tableSize);
		reserve(beta_b, tableSize);
		populate_alphabeta(hemi_f, beta_f, 1);
		populate_alphabeta(hemi_b, beta_b, 1);

		/*average backscattering attenuation*/
		reserve(A_b, tableSize);
		populate_A_b(A_b);
		
		
		/*average longitudinal shift*/
		reserve(delta_b, tableSize);
		populate_delta_b(delta_b);

		
		/*average backscattering deviation*/
		reserve(sigma_b, tableSize);
		populate_sigma_b(sigma_b);		
		
	}
	
    public void surface(output color Ci, Oi;)
    {
		// Get unit vectors along local axis in globle coordinate system
		vector lx  =  normalize(dPdu);
		vector ly  =  normalize(   N);  //the shading normal
		vector lz  =  normalize(dPdv);	//hair tangent (from root to tip)
		
        //I is the incident ray dir(from eye to the shading point)
		vector   omega_o = GlobalToLocal( -normalize(I), lx, ly, lz ); 
		float      phi_o = atan(omega_o[1], omega_o[0]);
		float    theta_o = M_PI_2 - acos(omega_o[2]);
		float alpha_back = radians(BackScattering_LongituShift);
		float  beta_back = radians(BackScattering_LongituWidth);
		
		illuminance(P) //P is the shading point position, a function of (u,v)
		{
			//L is light ray (from shading point to the light source)
			vector omega_i = GlobalToLocal( normalize(L), lx, ly, lz ); 
			float    phi_i = atan(omega_i[1], omega_i[0]);
			float    phi   = abs(phi_o - phi_i); //relative azimuth (within the normal plane)

			if ( phi > PI )
				phi -= 2 * PI;
			phi = abs(phi);

			float theta_i = M_PI_2 - acos(omega_i[2]);
			float theta   = theta_i + theta_o;
				  theta_h = theta * 0.5; //half longitudial angle (wrt the normal plane)

			//interpolate from pre-computation
			float Ab;
			float deltab;
			float sigmab;
			
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
			color f_direct_back  =  2 * Ab * g( sigmab + beta_back, theta_h - deltab + alpha_back) 
									/ (PI * cos(theta) * cos(theta));
				  f_direct_back = BackScattering_Color * BackScattering_Intensity * f_direct_back;

			color f_scatter_back = 2 * Ab * g( sigmab + sigma_f + beta_back, theta_h - deltab + alpha_back) 
									/ (PI * cos(theta) * cos(theta));
		          f_scatter_back = BackScattering_Color * BackScattering_Intensity * f_scatter_back;

			//single scattering for direct and indirect lighting
			color f_direct_s  =  M(  R, theta_h) * N_(  R, phi) 
							  +  M( TT, theta_h) * N_( TT, phi) 
							  +  M(TRT, theta_h) * N_(TRT, phi);
							  
			color f_scatter_s = MG(  R, theta_h) * NG(  R) 
							  + MG( TT, theta_h) * NG( TT) 
							  + MG(TRT, theta_h) * NG(TRT);
				  f_scatter_s = ForwardScattering_Color * ForwardScattering_Intensity * f_scatter_s;
			
			
			color F_direct  = illuminated * ( f_direct_s + f_direct_back );
			color F_scatter = (T_f - illuminated) * d_f * ( f_scatter_s + PI * d_b * f_scatter_back);

			//combine the direct and indirect scattering components
			Ci  += (F_direct + F_scatter) * cos(theta_i);
		}

		Oi  = Os; //Os is the surface opacity
		Ci *= Oi;
	
    }
}