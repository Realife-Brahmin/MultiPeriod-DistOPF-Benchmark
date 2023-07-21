function strArea = convert2doubleDigits(Area)
    switch Area < 10 %single digit number
        case true
            strArea = strcat("0", num2str(Area));
        case false
            strArea = num2str(Area);
        otherwise
            strArea = num2str(Area);
    end
end