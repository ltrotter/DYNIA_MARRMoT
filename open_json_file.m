function [thetas, OF, fluxes, stores] = open_json_file(file_name)

jsontext = fileread(file_name);

data = jsondecode(jsontext);

thetas = [data.theta];
OF = [data.OF];
fluxes = [data.fluxes];
stores = [data.stores];