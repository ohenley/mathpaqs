with Ada.Text_IO;

procedure Show_floats_limits is

  --  Presumably the Intel extended double precision (80 bits, i.e. 10 bytes)
  type Digits_18 is digits 18;

  generic
    type F is digits <>;
    name: String;
  procedure Show_limits;

  procedure Show_limits is
    use Ada.Text_IO;
  begin
    Put_Line(name);
    Put_Line("    Largest positive number  . . . . " & F'Image(F'Last));
    Put_Line("    Largest negative number  . . . . " & F'Image(F'First));
    Put_Line("    Significant decimal digits . . . " & Integer'Image(F'Digits));
    New_Line;
  end Show_limits;

  --  In many cases Float is the IEEE single precision type (32 bits, i.e. 4 bytes).
  --  Single precision accumulates quickly significant rouding errors.
  --  It should to be used only in special cases where memory is scarse.
  --
  procedure SF   is new Show_limits(Float,           "Float");

  --  Often, Long_Float is the IEEE dougle precision type (64 bits, i.e. 8 bytes).
  --
  procedure SLF  is new Show_limits(Long_Float,      "Long_Float");

  procedure SLLF is new Show_limits(Long_Long_Float, "Long_Long_Float");
  procedure SD18 is new Show_limits(Digits_18,       "Digits_18");

begin
  SF;
  SLF;
  SLLF;
  SD18;
end Show_floats_limits;
