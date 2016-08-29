function out = casadi2ampl(fname,nlp,data,discrete)
  import casadi.*
  fun = Function('nlp',nlp,char('x','p'),char('f','g'));
  fun = fun.expand();

  ins = fun.sx_in();
  fun = Function('f',ins,fun.call(ins),struct('live_variables',false));

  lbg = -inf*DM(fun.sparsity_out(1));
  ubg = inf*DM(fun.sparsity_out(1));
  if isfield(data,'lbg'), lbg(:,:) = data.lbg; end
  if isfield(data,'ubg'), ubg(:,:) = data.ubg; end


  lbx = -inf*DM(fun.sparsity_in(0));
  x0 = DM(fun.sparsity_in(0));
  ubx = inf*DM(fun.sparsity_in(0));
  if isfield(data,'lbx'), lbx(:,:) = data.lbx; end
  if isfield(data,'ubx'), ubx(:,:) = data.ubx; end
  if isfield(data,'x0'), x0(:,:) = data.x0; end
  
  p = DM(fun.sparsity_in(1));
  if isfield(data,'p'), p(:,:) = data.p; end
  
  main = fopen('results.mod','w');
  
  main = '';
  for i=1:size(p,1)
     main = [ main sprintf('param p%d = %.16f;\n',i-1,full(p(i)))];
  end
  for i=1:size(x0,1)
    if discrete(i)
        common = sprintf('var x%d := %.16f, integer',i-1,full(x0(i)));
    else
        common = sprintf('var x%d := %.16f',i-1,full(x0(i)));
    end
    if isinf(full(lbx(i))) && isinf(full(ubx(i)))
      main = [ main sprintf('%s;\n',common)];
    elseif isinf(full(lbx(i)))
      main = [ main sprintf('%s, <= %.16f;\n',common,full(ubx(i)))];
    elseif isinf(full(ubx(i)))
      main = [ main sprintf('%s, >= %.16f;\n', common,full(lbx(i)))];
    else
      main = [ main sprintf('%s, >= %.16f, <= %.16f;\n', common,full(lbx(i)),full(ubx(i)))];
    end
  end

  algorithm = strsplit(fun.getDescription(),'\n');
  for i=1:numel(algorithm)
    l = algorithm{i};
    if strfind(l, '@')
      break
    end
  end
  algorithm = strjoin(algorithm(i:end),'\n');
    
  algorithm = regexprep(algorithm,'@','at');
  algorithm = regexprep(algorithm,'input\[0\]\[(\d+)\]','x$1');
  algorithm = regexprep(algorithm,'input\[1\]\[(\d+)\]','p$1');
  algorithm = regexprep(algorithm,'output\[0\]\[(\d+)\]','f');
  algorithm = regexprep(algorithm,'output\[1\]\[(\d+)\]','g$1');
  algorithm = regexprep(algorithm,'sq\((.*?)\)','($1)^2');
  
  is_var = struct();
  
  algorithm_split = strsplit(algorithm, ';');
  
  for i=1:numel(algorithm_split)-1
    a = strtrim(algorithm_split{i});
    
    lhs_rhs = strsplit(a, '=');
    lhs = lhs_rhs{1};
    rhs = lhs_rhs{2};


    if ~isempty(strfind(rhs, '>')) || ~isempty(strfind(rhs, '<')) || ~isempty(strfind(rhs, '='))
      rhs = sprintf('if %s then 1 else 0', rhs);
    end
    var = ~isempty(strfind(rhs, 'x'));

    if ~var
      for m=regexp(rhs,'(at\d+)','match')
        m=m{1};
        if is_var.(m)
          var = true;
        end
      end
    end
    is_var.(strtrim(lhs)) = var;

    a = sprintf('%s = %s',lhs, rhs);
    if var
      main = [ main sprintf('var %s;\n', a)];
    else
      main = [ main sprintf('param %s;\n', a)];
    end
  end

  constr = '';
  for i=1:size(lbg,1)
    if isinf(full(lbg(i))) && isinf(full(ubg(i)))
      continue;
    elseif full(lbg(i))==full(ubg(i))
      constr = [constr sprintf('  con%d: g%i = %.16f;\n',i-1,i-1,full(lbg(i)))];
    elseif isinf(full(lbg(i)))
      constr = [constr sprintf('  con%d: g%i <= %.16f;\n',i-1,i-1,full(ubg(i)))];
    elseif isinf(full(ubg(i)))
      constr = [constr sprintf('  con%d: g%i >= %.16f;\n',i-1,i-1,full(lbg(i)))];
    else
      constr = [constr sprintf('  con%d_lower: g%i >= %.16f;\n',i-1,i-1,full(lbg(i)))];
      constr = [constr sprintf('  con%d_upper: g%i >= %.16f;\n',i-1,i-1,full(ubg(i)))];
    end
  end
  
  % Write model to file
  file = fopen([fname '.mod'],'w');
  fprintf(file,'reset;\n%s\nminimize f1: f;\ns.t.\n%s\n', main, constr);
  
  displayvar = cell(1,size(x0,1));
  for i=1:size(x0,1)
      displayvar{i} = sprintf('x%d',i-1);
  end
  
  
  file = fopen([fname '.com'],'w');
  fprintf(file,'solve;\ndisplay %s;\n',strjoin(displayvar,','));
    
end
