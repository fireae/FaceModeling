function R = Ang2Rot(angle)

    angle = angle.*(pi/180);
    rou = angle(1);
    phi = angle(2);
    psi = angle(3);
    
    R(1,1) = cos(psi)*cos(phi);
    R(1,2) = cos(psi)*sin(phi)*sin(rou) - sin(psi)*cos(rou);
	R(1,3) = cos(psi)*sin(phi)*cos(rou) + sin(psi)*sin(rou);

	R(2,1) = sin(psi)*cos(phi);
	R(2,2) = sin(psi)*sin(phi)*sin(rou) + cos(psi)*cos(rou);
	R(2,3) = sin(psi)*sin(phi)*cos(rou) - cos(psi)*sin(rou);

	R(3,1) = -sin(phi);
	R(3,2) = cos(phi)*sin(rou);
	R(3,3) = cos(phi)*cos(rou);
end  