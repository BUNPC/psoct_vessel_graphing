function parsave(fname, var)
    orignal_name = inputname(2);
    S.(orignal_name) = var;
    save(fname,  '-struct', 'S')
end