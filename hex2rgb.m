function rgb = hex2rgb(hex)
    hex = hex(2:end); % Remove '#'
    rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))] / 255;
end