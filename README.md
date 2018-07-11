# Solar-Shock-Diagnosis

Background introduction to the shock in the solar atmosphere:

The Sun is always under oscillation on its surface. The ever-present osciallation is driven by its internal activities including thermalnuclear fusion at its core, collisional impact in the convection zone, and so on.
The oscillation launches compressible waves into the solar atmosphere above the surface (photosphere). As propagating in the straitified solar atmosphere, the compressible waves become large amplitude and nonlinearly evolve into shocks. The shocks are then dissipated and play a crucial role in heating the solar atmosphere. Therefore, it is an important issue to comprehensively diagnose the shocks in the solar atmosphere.


Advanced observations of the shock activities from IRIS:

IRIS, a short name standing for "Interface Region Imaging Spectrograph", is a latest space-based solar observatory. It is featured by its pioneering optical design, which allows for simultaneous 2D imaging and spectroscopic observation by placing a mirror with a slit at the focus position. Two optical sub-systems extend after the first imaging location: one is the re-imaging system to record the image as refelected from the mirror, and the other is the spectrograph system to record the emission spectral profiles passing through the slit. How to utilize such observations to diagnose the shock status is the goal of this project.


An approach to comprehensively diagnose solar atmosphere shocks:

This method self-consistently incorporates the theoretical Rankine-Hugoniot jump conditions and the quantities derived from simultaneous imaging and spectroscopic observations. According to the algebraic derivation, we have a set of four equations with the following information as the known conditions: emission intensities (I1, I2) and Doppler shifts (D1, D2) upstream and downstream of the shock, projected propagation speed in the plane of sky (V_beta), temperature response function for a certain band of the spectrograph (G(T)). This set of equations can ultimately be transformed into three expressions of M1 (upstream Mach number), T1 (upstream temperature), T2 (downstream temperature) with V_alpha (propagation speed in the line-of-sight direction) and an equation with only V_alpha to be solved iteratively. 

By employing this method to the observations from Interface Region Imaging Spectrograph (IRIS), the full set of shock wave's parameters can be derived, such as the bulk velocity and temperature of the plasma in the upstream and downstream, the propagation speed and direction. This method provides an effective tool in comprehensive diagnosing shock waves and thereby studying the heating and the temperature profiles of solar atmosphere.


How to use the code in this project to the diagnosis of solar atmosphere shock:

Please refer to the document "???.pdf" for the details, and run the code "???.py" for an example test. The input parameters in the code "???.py" can be adjusted according to the user's case. As the output, the full set of shock paramters can be calculated from the code "???.py".

Hope this project help you in studying the solar atmosphere shocks.


Jiansen HE
Group leader of the project
jshept@pku.edu.cn
Peking University

July 11, 2018

