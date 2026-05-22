// THIS IS A QUESTION FOR THE TA'S, WHY ISNT THIS EXACTLY DELTA_E???
double get_energy_n(double spin[N][N][3]){
    double Energy_n = 0;

    for (int di = -NC; di <= NC; di++) {
        int i = (n1 + di + N) % N;
        for (int dj = -NC; dj <= NC; dj++) {
            int j = (n2 + dj + N) % N;

            // Magnetic Field energy
            double E_H = -Hz * spin[i][j][2];

            // indexes for spin spin interactions
            int i_2 = (i + 1 + N) % N;
            int i_3 = (i - 1 + N) % N;
            int j_2 = (j + 1 + N) % N;
            int j_3 = (j - 1 + N) % N;


            // spin spin alinging interaction
            double E_J = 0;
            E_J += -J/2*dot_prod(spin[i][j],spin[i_2][j]);
            E_J += -J/2*dot_prod(spin[i][j],spin[i_3][j]);
            E_J += -J/2*dot_prod(spin[i][j],spin[i][j_2]);
            E_J += -J/2*dot_prod(spin[i][j],spin[i][j_3]);

            // spin orbit coupling
            double E_D = 0;
            E_D += -D/2*cross_x(spin[i][j],spin[i_2][j]);
            E_D += -D/2*cross_x(spin[i][j],spin[i_3][j]);
            E_D += -D/2*cross_y(spin[i][j],spin[i][j_2]);
            E_D += -D/2*cross_y(spin[i][j],spin[i][j_3]);

            Energy_n += E_H + E_J + E_D;
        }
    }

    return Energy_n;
}
