var radius = 50;

    macro "Circle Tool - C00cO11cc" {
        getCursorLoc(x, y, z, flags);
        makeOval(x-radius, y-radius, radius*2, radius*2);
    }

    macro "Circle Tool Options" {
        radius = getNumber("Radius: ", radius);
    }