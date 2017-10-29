function [I1,I2] = generate_test_images(n)

testcase = 2;

[X,Y] = meshgrid(1:n,1:n);

switch testcase
    case 1
        %% one Gaussian moving to the lower right
        position = [0.48,0.49];
        sigma = 0.15;
        shift = [0.04,0.02];
        
        abs_pos_y = position(1)*n;
        abs_pos_x = position(2)*n;
        sigma_scaled = sigma*n;
        
        I1 = 255*exp(-(((Y-abs_pos_y).^2 + (X-abs_pos_x).^2))/(2*sigma_scaled^2));
        
        abs_pos_y_fin = (position(1)+shift(1))*n;
        abs_pos_x_fin = (position(2)+shift(2))*n;
        
        I2 = 255*exp(-(((Y-abs_pos_y_fin).^2 + (X-abs_pos_x_fin).^2))/(2*sigma_scaled^2));
    case 2
        %% two Gaussians circling around the center
        
        % construction of first Gaussian
        position = [0.3,0.5];
        sigma = 0.05;
        shift = [0.05,0.05];
        
        abs_pos_y = position(1)*n;
        abs_pos_x = position(2)*n;
        sigma_scaled = sigma*n;
        
        I1_1 = 255*exp(-(((Y-abs_pos_y).^2 + (X-abs_pos_x).^2))/(2*sigma_scaled^2));

        abs_pos_y_fin = (position(1)+shift(1))*n;
        abs_pos_x_fin = (position(2)+shift(2))*n;
        
        I2_1 = 255*exp(-(((Y-abs_pos_y_fin).^2 + (X-abs_pos_x_fin).^2))/(2*sigma_scaled^2));
        
        % construction of second Gaussian
        position = [0.7,0.5];
        sigma = 0.10;
        shift = [-0.05,-0.05];
        
        abs_pos_y = position(1)*n;
        abs_pos_x = position(2)*n;
        sigma_scaled = sigma*n;
        
        I1_2 = 255*exp(-(((Y-abs_pos_y).^2 + (X-abs_pos_x).^2))/(2*sigma_scaled^2));

        abs_pos_y_fin = (position(1)+shift(1))*n;
        abs_pos_x_fin = (position(2)+shift(2))*n;
        
        I2_2 = 255*exp(-(((Y-abs_pos_y_fin).^2 + (X-abs_pos_x_fin).^2))/(2*sigma_scaled^2));
        
        % take the maximum of the two Gaussians
        I1 = max(I1_1,I1_2);
        I2 = max(I2_1,I2_2);
end