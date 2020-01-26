function fitnessvalue = fitnessf(x)
%HILLFIT 此处显示有关此函数的摘要
%   此处显示详细说明
%    y1 = -15:0.1:15;
%    y2 = -15:0.1:15;
global FUN_NUM FITCOUNT Fvalues
[m, D] = size(x);
FITCOUNT = FITCOUNT + m;
switch FUN_NUM
    case 1
        fitnessvalue=fitnessOLD( x,1);
    case 2
        fitnessvalue=FitnessCBG( x );
    case 3
        fitnessvalue=1-1*FitnessMSG( x )';
    case 4
        fitnessvalue=FitnessTFG( x, 1 );
    case 5
        fitnessvalue=FitnessNGLI( x ); 
    case 6
        fitnessvalue=FitnessCLUST( x ); 
    case 7
        fitnessvalue=FitnessTTG( x ); 
    case 8
        fitnessvalue=FitnessRFC( x ); 
    case 9
        fitnessvalue=FitnessNPEAKS( x ); 
    case 10
        fitnessvalue=FitnessNM( x ); 
        
    case 11
        fitnessvalue = fgeneric(x')' - fgeneric('ftarget');
end

Fvalues = [Fvalues;fitnessvalue];

end

