%% panelcode(airfoil, angle of attack (º), airflow velocity (m/s), resolution)

figure()
panelcode("naca633618.dat",10,2,50) % High Res
figure()
panelcode("naca633618.dat",-10,2,50) % Downward aoa, displaying pressure variation 
figure()
panelcode("n0012.dat", 0,2,50) % Symmetrical airfoil
