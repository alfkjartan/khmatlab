function sgm=switchalias(sgml,al)
% Substitutes the values in sgml with the aliases in al

% Kjartan Halvorsen
% 2001-03-05

[rows,cols]=size(sgml);

sgm=cell(rows,cols);

for c=1:cols
   for r=1:rows
      sgm{r,c}=getvalue(al,sgml{r,c});
   end
end
