function angle = Rot2Ang(R)

angle(1) = atan2(R(3,2),R(3,3));
angle(2) = atan2(-R(3,1),sqrt(R(3,2)^2 + R(3,3)^2));
angle(3) = atan2(R(2,1),R(1,1));
end