package body Discrete_Random_Simulation is

  function Index_Linear_Search (
    U01 : Probability_value;  --  Probability value. For simulation: random, uniform in [0,1]
    Fx  : Probability_array   --  Fx is the Cumulative distribution function (CDF), F(x)
  )
  return Integer
  is
  begin
    --  Search from top to bottom:
    for i in reverse Fx'Range loop
      if Fx(i) <= U01 then
        return i;
      end if;
    end loop;
    return Fx'First;
  end Index_Linear_Search;

  function Index_Dichotomic_Search (
    U01 : Probability_value;  --  Probability value. For simulation: random, uniform in [0,1]
    Fx  : Probability_array   --  Fx is the Cumulative distribution function (CDF), F(x)
  )
  return Integer
  is
    l, r, i: Integer;
  begin
    l:= Fx'First;
    r:= Fx'Last;
    loop
      i:= l + ((r - l) / 2);  --  This is (l+r) / 2 without overflow problem
      if Fx(i) <= U01 then
        l:= i;
      else
        r:= i;
      end if;
      exit when r-l <= 1;
    end loop;
    --  Finally, a top-bottom search on Fx(r) >= Fx(i) >= Fx(l)
    if Fx(r) <= U01 then
      i:= r;
    elsif Fx(i) <= U01 then
      null;
    else
      i:= l;
    end if;
    return i;
  end Index_Dichotomic_Search;

  function To_cumulative (p: Probability_array; check: Boolean:= False) return Probability_array is
    sum: Probability_value := 0.0;
    Fx: Probability_array(p'Range);
    tol: constant := 1.0e-7;
  begin
    for i in p'Range loop
      Fx(i):= sum;  --  CDF begins with 0.
      sum := sum + p(i);
      if check and then p(i) not in 0.0 .. 1.0 then
        --  This check is redundant if the type Real is aleady a subtype in 0.0 .. 1.0
        raise Constraint_Error with "Probability value not in [0,1]";
      end if;
    end loop;
    if check and then abs(sum-1.0) > tol then
      raise Constraint_Error
        with "Sum of probabilities should be 1, is" & Probability_value'Image(sum);
    end if;
    return Fx;
  end To_cumulative;

end Discrete_Random_Simulation;
