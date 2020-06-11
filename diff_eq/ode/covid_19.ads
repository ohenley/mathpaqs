package COVID_19 is

   type Real is digits 15;

   type Scenario is (No_Lockdown, Lockdown, Lockdown_in_two_Steps);

   type Status is (Susceptible, Exposed, Infectious, Recovered);

   type Simulation_Data is array (Status range <>, Natural range <>) of Real;

   function Simulation (s: Scenario; n_population: Integer; n_iter : Integer) return Simulation_Data;

   procedure Export_To_CSV (sim_data: Simulation_Data; s: Scenario; n_population: Integer; n_iter : Integer);

end COVID_19;
