
for(i=0;i<t_cs;i++){
		// Calculate reaction rates and enthalpy changes
		Q_t=0;
		/*for(j=0;j<n_rs;j++){
			k_f=1; k_r=1;
			for(l=0;l<n_cs;l++){
				k_f*=pow(psi[l*t_cs+i],stoi[j*n_cs+l]);
				k_r*=pow(psi[l*t_cs+i],stoi[n_cs*n_rs+j*n_cs+l]);
			}
			R_f[j]=A_f[j]*exp(-E_f[j]/eps[i])*omega[k]*k_f;
			R_r[j]=A_r[j]*exp(-E_r[j]/eps[i])*omega[k]*k_r;
			Q_t+=del_H_r[j]*(R_f[j]-R_r[j]);
		}*/
		
		// Collision operations and concentration changes due to reactions
		for(k=0;k<9;k++){
			// Fluid
			f[k*t_cs+i]+=f_0[k*t_cs+i]*(i_tau_v+b_g0*i_T_0*(eps[i]-T_0)*c*(e_is[2*k+1]-u[t_cs+i]))-i_tau_v*f[k*t_cs+i];
			//f[k*t_cs+i]+=i_tau_v*(f_0[k*t_cs+i]-f[k*t_cs+i]);		// Convection off
            
			// Internal energy
			f[(k+9)*t_cs+i]+=i_tau_c*(f_0[(k+9)*t_cs+i]-f[(k+9)*t_cs+i])+Q_t;							// Collision and heats of reaction
            
			// Chemical species
			for(l=0;l<n_cs;l++){
				f[((2+l)*9+k)*t_cs+i]+=i_tau_s[l]*(f_0[((2+l)*9+k)*t_cs+i]-f[((2+l)*9+k)*t_cs+i]);	// Collision
				/*for(j=0;j<n_rs;j++){
					f[((2+l)*9+k)*t_cs+i]-=stoi[j*n_cs+l]*R_f[j]-stoi[n_cs*(n_rs+j)+l]*R_r[j];		// Reactions
				}*/
			}
		}
	}
}