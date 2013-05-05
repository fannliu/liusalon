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
	constant float M_1_PI = 1.0 / PI;
	constant float M_2_PI = 2.0 / PI;
	
	constant float segment = 0.1;
	constant float inv_segment = 10;
	
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
	color M_(uniform float component; float theta_h;)
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
	
	
	color fs(float theta, phi;)
	{
		float theta_h = theta * 0.5;
		color f_R   =   PrimaryHL_Color *   PrimaryHL_Intensity * M_(R, theta_h)   * N_(R, phi);
		color f_TT  =  BacklitRim_Color *  BacklitRim_Intensity * M_(TT, theta_h)  * N_(TT, phi);
		color f_TRT = SecondaryHL_Color * SecondaryHL_Intensity * M_(TRT, theta_h) * N_(TRT, phi);

		return (f_R + f_TT + f_TRT)/(cos(theta) * cos(theta));
	}
	
	color integrateOverFullSphere()//integrate over the full sphere around the shading point
	{
		uniform float theta, phi;
		uniform float result = 0.0;
		
		for(phi = - PI; phi <= PI; phi += segment){
			for(theta = - PI; theta <= PI; theta += segment)
				result  += fs(theta, phi);
				
		
		result *= segment;
		
		return result;
		
	}
	
	color fs_normalized(float theta, phi;)
	{
		return fs(theta, phi)/integrateOverFullSphere();
		
	}
	color color_abs(uniform color x){
		return color(abs(x[0]), abs(x[1]), abs(x[2]));
	}
	
	color color_pow(uniform color x; uniform float n;){
		return color(pow(x[0],n), pow(x[1],n), pow(x[2],n));
	}
	
	color interpolate_theta_i(varying color ary[]; varying float theta_i) {
        if (theta_i == - M_PI_2)
            return ary[0];
        else if (theta_i < M_PI_2)
            return ary[0] - (ary[1] - ary[0]) * (M_PI_2- theta_i) * inv_segment;
        else if (theta_i == M_PI_2)
            return ary[tableSize-1];
        else if (theta_i > M_PI_2)
            return ary[tableSize-1] + (ary[tableSize-1] - ary[tableSize-2]) * (theta_i - M_PI_2) * inv_segment;
        else {
            float offset = (theta_i - M_PI_2) * inv_segment;
            float low = floor(offset);
            float high = ceil(offset);
            if (low == high)
                return ary[low];
            else
                return ary[low] + (ary[high] - ary[low]) * (offset - low)/(high - low);
        }
    }
	
	float isHemisphere(uniform float hemisphere, theta_i, phi_i, theta_o, phi_o;)
    {
        // map phi_i from [-PI, PI] to [0, 2*PI]
        if (phi_i >= 0)
            phi_i -= M_PI_2;
        else
            phi_i += 3*M_PI_2;
			
        uniform float cos_theta_i = cos(theta_i);
        uniform vector vi = (-sin(phi_i)*cos_theta_i, cos(phi_i)*cos_theta_i, sin(theta_i));
 
        // map phi_o from [-PI, PI] to [0, 2*PI]
        if (phi_o >= 0)
            phi_o -= M_PI_2;
        else
            phi_o += 3*M_PI_2;
        uniform float cos_theta_o = cos(theta_o);
        uniform vector vo = (-sin(phi_o)*cos_theta_o, cos(phi_o)*cos_theta_o, sin(theta_o));
 
        uniform float IdotO = vi . vo;
        // need front hemisphere
        if (hemisphere == hemi_f) { 
            if (IdotO < 0)
                return 0;
            else
                return 1;
        }
        else {//hemi_b
            if (IdotO >= 0)
                return 0;
            else
                return 1;
        }
    }
    
	//pre-computing and tabulating the uniform variables in the Constructor
	void populate_a(uniform float hemisphere; output uniform color a[]) 
    {
		uniform float theta_i = - M_PI_2;
		uniform float phi_o, theta_o, phi_i;
		
        uniform float i;
        
        for (i = 0; i < tableSize; i += 1) {
 
            uniform color sum_phi_o = 0.0;
			
			//omega_o decompose to phi_o and theta_O
            for (phi_o = - PI; phi_o <= PI; phi_o += segment) {
 
                uniform color sum_theta_o = 0.0;
				
                for (theta_o = - M_PI_2; theta_o <= M_PI_2; theta_o += segment) {
 
					float theta = theta_i + theta_o;
					
					uniform color sum_phi_i = 0.0;
                    
					for (phi_i = - PI; phi_i <= PI; phi_i += segment) {
                        if (isHemisphere(hemisphere, theta_i, phi_i, theta_o, phi_o) == 0) 
                            continue;
                        
						float phi = abs(phi_o - phi_i);
						if( phi > PI)
							phi -= 2*PI;
						phi = abs(phi);
						
						
                        uniform color fcos = color_abs(fs_normalized(theta, phi)) * cos(theta_i);
                        sum_phi_i += fcos;
                    }
					
                    sum_phi_i *= segment;
                    sum_theta_o += sum_phi_i;
                }
                sum_theta_o *= segment;
                sum_phi_o += sum_theta_o;
            }
            sum_phi_o *= segment;
			
            push( a, sum_phi_o * M_1_PI );
            theta_i += segment;
        }
    }
	
	color integrateOverHemisphere(uniform float theta_i; uniform float hemisphere)
	{
	    
		uniform float phi_o, theta_o, phi_i;

        uniform color sum_phi_o = 0.0;
		
		uniform float inv_fs_integral = 1.0 /integrateOverFullSphere();
		
       for (phi_o = - PI; phi_o <= PI; phi_o += segment) {
 
            uniform color sum_theta_o = 0.0;
			
            for (theta_o = - M_PI_2; theta_o <= M_PI_2; theta_o += segment) {
                
				float theta = theta_i + theta_o;
				uniform color sum_phi_i = 0.0;
					
                for (phi_i = - PI; phi_i <= PI; phi_i += segment) {
                    if (isHemisphere(hemisphere, theta_i, phi_i, theta_o, phi_o) == 0) 
                        continue;
						
                    float phi = abs(phi_o - phi_i);
					if( phi > PI)
						phi -= 2*PI;
					phi = abs(phi);	

                    uniform color f = fs_normalized(theta, phi);
					
                    sum_phi_i  += f;
                }
                sum_phi_i  *= segment;
                sum_theta_o += sum_phi_i;
            }
            sum_theta_o *= segment;
            sum_phi_o += sum_theta_o;
        }
        sum_phi_o *= segment;
 
        return sum_phi_o;
	}
	
	color integrateOverHemisphereWeighted(uniform float theta_i; uniform float coef; uniform float hemisphere){

		
		uniform float coef_R, coef_TT, coef_TRT;
        if (coef == 0) {
            coef_R = radians(PrimaryHL_LongituShift);
            coef_TT = radians(BacklitRim_LongituShift);
            coef_TRT= radians(SecondaryHL_LongituShift);
        }
        else {
            coef_R = radians(PrimaryHL_LongituWidth);
            coef_TT = radians(BacklitRim_LongituWidth);
            coef_TRT = radians(SecondaryHL_LongituWidth);
        }
 
        uniform float phi_o, theta_o, phi_i;
		
		uniform float inv_fs_integral = 1.0 /integrateOverFullSphere();
		
		for (phi_o = - PI; phi_o <= PI; phi_o += segment) {
 
            uniform color sum_theta_o = 0.0;
			
            for (theta_o = - M_PI_2; theta_o <= M_PI_2; theta_o += segment) {
                
				float theta = theta_i + theta_o;
				float theta_h = theta * 0.5;
				float inv_cos_theta_pow2 = 1.0 / pow(cos(theta),2);
				
				uniform color sum_phi_i = 0.0;
					
                for (phi_i = - PI; phi_i <= PI; phi_i += segment) {
                    if (isHemisphere(hemisphere, theta_i, phi_i, theta_o, phi_o) == 0) 
                        continue;
						
                    float phi = abs(phi_o - phi_i);
					if( phi > PI)
						phi -= 2*PI;
					phi = abs(phi);	

                    color fs_R   =   PrimaryHL_Color *   PrimaryHL_Intensity * M_(R, theta_h)   * N_(R, phi);
					color fs_TT  =  BacklitRim_Color *  BacklitRim_Intensity * M_(TT, theta_h)  * N_(TT, phi);
					color fs_TRT = SecondaryHL_Color * SecondaryHL_Intensity * M_(TRT, theta_h) * N_(TRT, phi);

                    uniform color f = (fs_R * coef_R + fs_TT * coef_TT + fs_TRT * coef_TRT) * inv_cos_theta_pow2 * inv_fs_integral;
					
                    sum_phi_i  += f;
                }
                sum_phi_i  *= segment;
                sum_theta_o += sum_phi_i;
            }
            sum_theta_o *= segment;
            sum_phi_o += sum_theta_o;
        }
        sum_phi_o *= segment;
 
        return sum_phi_o;
		
	}
	
	void populate_alphabeta(uniform float hemisphere; output uniform color target[]; float coef;)
    {
        uniform float theta_i = - M_PI_2;
        uniform float i;
        uniform color denominator = 0.0;
        uniform color numerator = 0.0;
 
        for (i = 0; i < tableSize; i++) {
			  numerator = integrateOverHemisphereWeighted(theta_i, coef, hemisphere);
            denominator = integrateOverHemisphere(theta, hemisphere);
             
            push( target, numerator / denominator );
            theta_i += segment;
        }
    }

	void populate_A_b(output uniform color A[])
    {
        uniform float i;
 
        for (i = 0; i < tableSize; i += 1) {
            uniform color af = a_f[i];
            uniform color ab = a_b[i];
			
			uniform color afPow2 = color_pow(af, 2);
			
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
			
			uniform color afPow2 = color_pow(af, 2);
			uniform color abPow2 = color_pow(ab, 2);
			
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
			
			uniform color betabPow2 = color_pow(betab, 2);
			uniform color betafPow2 = color_pow(betaf, 2);
			uniform color abPow3 = color_pow(ab,3);
			
			uniform color sigmab_r = (1 + d_b * pow(af,2)) * (ab + abPow3) * sqrt(2*betafPow2 + betabPow2)
									/ ( ab + abPow3*(2*betaf + 3*betab) );
			/*uniform float sigmab_g = (1 + d_b * pow(af[1],2)) * (ab[1] + abPow3[1]) * sqrt(2*betafPow2[1] + betabPow2[1])
									/ ( ab[1] + abPow3[1]*(2*betaf[1] + 3*betab[1]) );	*/
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

        // get deep shadow map value
        float shadow_bias = 0.005;
        float shadow_blur = 0.002;
        float shadow_samples = 9;
        uniform string shadowmap_path = "";
        attribute("light:user:delight_shadowmap_name", shadowmap_path);

        float shadow_p = shadow(shadowmap_path, P, "bias", shadow_bias, "samples", shadow_samples, "blur", shadow_blur);
 		
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

			//interpolate value from pre-computation
			float     Ab = interpolate_theta_i(    A_b, theta_i);
			float deltab = interpolate_theta_i(delta_b, theta_i);
			float sigmab = interpolate_theta_i(sigma_b, theta_i);
			
			// compute the amount of shadow from the deep shadow maps
            float shadowed = shadow_p;
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
			color f_direct_s  =  M_(  R, theta_h) * N_(  R, phi) 
							  +  M_( TT, theta_h) * N_( TT, phi) 
							  +  M_(TRT, theta_h) * N_(TRT, phi);
							  
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