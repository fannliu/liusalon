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
	//Defining the uniform parameters unsed
	float R   = 0;
	float TT  = 1;
	float TRT = 2
	//Defining the varying parameters which have different values for each shading point
	varying float hairs_in_front;
	varying float phi;
	varying float theta_h;
	
    float g(float deviation, x;)
    {	//unit-integral zero-mean Gaussian distribution
       return exp( - x*x /( 2*deviation*deviation ) ) / ( deviation * sqrt(2*PI) );
    }

	float N(component){
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

	float NG(component){
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

	color M(component)
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

	color MG(component)
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

	color compute_f_direct_back(){
		color f_direct_back = 0;//need to fill
		return BackScattering_Color * BackScattering_Intensity * f_direct_back;
	}
	color compute_f_scatter_back(){
		color f_scatter_back = 0;//need to fill
		return BackScattering_Color * BackScattering_Intensity * f_scatter_back;
	}

	color compute_f_direct_s(){
		return M(R)*N(R) + M(TT)*N(TT) + M(TRT)*N(TRT);
	}

	color compute_f_scatter_s(){
		color f_scatter_s = MG(R)*NG(R) + MG(TT)*NG(TT) + MG(TRT)*NG(TRT);
		return ForwardScattering_Color * ForwardScattering_Intensity * f_scatter_s;
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
	//to-do: pre-computing and tabulating the uniform variables in the contructors
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

		
		illuminance(P) //P is the shading point position, a function of (u,v)
		{
			vector omega_i = GlobalToLocal( normalize(L), lx, ly, lz ); //L is light ray (from shading point to the light source)
			float phi_i    = atan(omega_i[1], omega_i[0]);
				  phi      = abs(phi_o - phi_i); //relative azimuth (within the normal plane)

			if ( phi > PI )
				phi -= 2 * PI;
			phi = abs(phi);

			float theta_i = PI * 0.5 - acos(omega_i[2]);
			float theta   = theta_i + theta_o;
				  theta_h = theta * 0.5; //half longitudial angle (wrt the normal plane)

			//compute the amount of shadow from the deep shadow maps?????


			color f_direct_back  = compute_f_direct_back();
			color f_scatter_back = compute_f_scatter_back();

			color f_direct_s  = compute_f_direct_s();
			color f_scatter_s = compute_f_scatter_s();
			
			
			color F_direct = directFraction * ( f_direct_s + f_direct_back);
			color F_scatter = (T_f - directFraction) * d_f * ( f_scatter_back + f_scatter_s);

			Ci  += (F_direct + F_scatter) * cos(theta_i);
		}

		Oi  = Os; //Os is the surface opacity
		Ci *= Oi;
	
    }
}