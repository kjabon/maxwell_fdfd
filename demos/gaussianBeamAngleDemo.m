%demo7 script
%how is "far field" determined? at what z?

filename = 'gaussianDemo.gif'
ns = .2:.05:2;
for n = ns
    close all
    h = figure;
    hold on;
    
    demo7(n)
    disp(n);
    
    lambda = 1;
    w0 = n*lambda; 
    z = 0:100;
    w = (w0^2*(1+(lambda*z/(pi*w0^2)).^2)).^(1/2);
    w1 = w+100+w0;
    w2 = 100-w0-w;
   
    y = 100:200;
    y2 = 100:-1:0;
    plot(w1,y);
    plot(w2,y);
    plot(w1,y2);
    plot(w2,y2);
    
    ylim([0,200])
    xlim([0,200])
    
%     theta = [1/(pi*n), -1/(pi*n)];
%     slope = 1./tan(theta);
%     y = [];
%     x = 1:200;
%     for j = 1:length(slope)
%         y(j,:) = slope(j)*x-100*slope(j)+100;
%         plot(x,y(j,:));
%     end
    drawnow;
    
    frame = getframe(h); 
    image = frame2im(frame);
    [indexedImg, colorMap] = rgb2ind(image, 256);

    if n == ns(1)

        imwrite(indexedImg, colorMap, filename, 'gif', 'Loopcount', inf,'DelayTime',0);
    else
        imwrite(indexedImg, colorMap, filename, 'gif', 'WriteMode', 'append','DelayTime',0);
    end
end
    