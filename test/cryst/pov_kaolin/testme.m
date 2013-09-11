#! /etc/alternatives/octave -q

addpath("../../../src");

cr = struct();
cr.name = "Kaolinite";
cr.a = [9.89398 9.86348 14.21688];
cr.b = [84.12595 98.99487 59.875919] * pi /180;
cr.nat = 17;
cr.ntyp = 4;
cr.x = [
        0.789999846   0.196788311   0.476933402 
        0.118263237   0.540299055   0.474753234 
        0.674705640   0.658679906   0.095081671 
        0.331910703   0.344666285   0.093717515 
        0.400339870   0.306418047   0.315410143 
        0.781894814   0.543270395   0.316464234 
        0.495496714   0.503977765   0.000839449 
        0.439287257   0.025507651   0.024902194 
        0.964007519   0.563917058   0.000936144 
        0.014581870   0.921329376   0.327708920 
        0.505627403   0.434897303   0.602411062 
        0.893580335   0.821285621   0.609049550 
        0.117848627   0.213566794   0.602543810 
        0.200245576   0.921187845   0.327094474 
        0.223705136   0.117829584   0.730325104 
        0.545975044   0.459137628   0.727808240 
        0.856830410   0.773541375   0.726620854 
        ];
cr.typ = [1 1 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4];
cr.attyp = {"Al","Si","O","H"};
cr.ztyp = [13 14 8 1];

mol = cr_crystalbox(cr,[-0.33 -0.33 -0.20],[1.30 1.30 1.80]);

## automatically add all the balls
# rep = mol_ball(mol)

# add the balls by hand
rep = mol_ball(mol,:,"Al",0,0.6,[127 127 192]);
rep = mol_ball(mol,rep,"Si",0,0.5,[89 89 89]);
rep = mol_ball(mol,rep,"O",0,0.3,[255 0 0]);
rep = mol_ball(mol,rep,"H",0,0.17,[230 230 230]);

## automatically add the sticks
#rep = mol_stick(mol,rep,:,:,[-1 1.05]);

# add the sticks by hand
rep = mol_stick(mol,rep,"Al","O",[1.5 2.0],0,0.02,[0 0 127]);
rep = mol_stick(mol,rep,"Si","O",[1.5 1.8],0,0.02,[0 0 127]);
rep = mol_stick(mol,rep,"O","H",[0.8 1.1],0,0.02,[0 0 127]);
rep = mol_stick(mol,rep,"O","H",[1.8 2.3],0,0.02,[0 127 0]);

# set camera, lights and write the pov
rep = rep_setdefaultscene(rep);
rep_write_pov(rep,"kaolin.pov");
