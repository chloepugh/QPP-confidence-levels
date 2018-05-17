;Model to fit a broken power law in log space

FUNCTION broke_pow_law,x,p

model = -(x LT p[1])*p[0]*x + $
         (x GT p[1])*(-(p[0]-p[3])*p[1] - p[3]*x) + p[2]

RETURN, model

END
