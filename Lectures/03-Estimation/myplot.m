function myplot

  bounds = [ -2, 8 ];

  y1 = [];
  y2 = [];
  x = [];
  for p = bounds(1):.1:bounds(2)
    x = [ x, p ];
    y1 = [ y1, f(p) ];
    y2 = [ y2, fp(p) ];
  end

  h0 = plot(bounds,[0, 0], 'k-');
  set(h0(1),"linewidth",5);
  hold on;
  h1 = plot(x,y1,'r-');
  set(h1(1),"linewidth",5);
  h2 = plot(x,y2,'g-');
  set(h2(1),"linewidth",5);
  title('Newton''s method for F(p) = (2p-4)(p^2 - 4p + 5), p0 = 4');
  axis on;
  hold off;

end

function y = f(p)
  y = ( 2 * p - 4 ) * ( p^2 - 4*p + 5 );
end

function y = j(p)
  y = (2 * p - 4)^2 + 2 * (p^2-4*p+5);
end

function y = fp(p)
  p0 = 4;
  y = f(p0) + j(p0) * (p - p0);
end

