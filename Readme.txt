I develop the PassiveFW package mainly for modeling active and passive seismic full wavefields with small scale geometry in layered media.
This package is mainly based on stiffness matrix method and discrete wavenumber method. 
main_active_shotgather.m show an example of active seismic wavefields in layered media due to a point load of Ricker wavelet.
main_passive_noise.m give an example of ambient noise with roadside geometry.


Just remember the time step dt is the key parameter that affect calculation speed.
A typical roadside example takes 50 minutes on a 4-core computer, but less than ten minutes on a 40-core computer. 

I will add more detailed instructions once I finish it.

For any suggestion, question or bug-fix, please don't hesitate to contact me.
Kai Zhang
Institute of Geophysical and Geochemical Exploration, 
Chinese Academy of Geological Sciences
84 Jinguang Road£¬Guangyang District
Hebei Langfang  065000, China	
Email:  naturekai@126.com


If you use this program in your research work, I suggest you cite the following article:

1.Kai Zhang, Hongyi Li, Xiaojiang Wang, and Kai Wang, Retrieval of shallow S-wave profiles from seismic reflection surveying and traffic-induced noise£¬Geophysics, 2020£¬85(6), EN105¨CEN117.

2.Bouchon, M., 2003, A Review of the Discrete Wavenumber Method: Pure and Applied Geophysics, 160, no. 3, 445-465. 

3.Kausel, E., 2018, Generalized stiffness matrix method for layered soils: Soil Dynamics and Earthquake Engineering, 115, 663-672.

