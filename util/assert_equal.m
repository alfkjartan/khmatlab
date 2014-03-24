function ok = assert_equal(str, a, b, thr)
  if ~max(max(abs(a-b)>thr))
    disp([str, '............OK'])
  else
    disp([str, '............FAILURE'])
    disp('Expected')
    a
    disp('Found')
    b
  end
end
