function result=test_iData_pca

  A = [269.8 38.9 50.5
272.4 39.5 50.0
270.0 38.9 50.5
272.0 39.3 50.2
269.8 38.9 50.5
269.8 38.9 50.5
268.2 38.6 50.2
268.2 38.6 50.8
267.0 38.2 51.1
267.8 38.4 51.0
273.6 39.6 50.0
271.2 39.1 50.4
269.8 38.9 50.5
270.0 38.9 50.5
270.0 38.9 50.5
];

  a=iData(A);
  b = pca(a);
  if all(abs(sum(b.Coefficients) - [ -0.6142    1.6051 ])  < 1e-4 )
    result = 1;
  else 
    result = 0; end
  
