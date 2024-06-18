classdef variables_de_inicio

    properties
      alfainicio = 20
      b = 10
      ds = 2
      dsta = 2
      fibdens = 1.77
      fmc = 60
      fvc = 0
      niu = 0.14
      NP = 6
      NR = 10
      resdens = 1.2
      rfm = 15
      rim = 15
      TEX = 800
      vel = 1.5
      zfg = 300
      zfm = 1000
      zi = 0
      zig = 0
      zim = -500
    end

    properties (Dependent)
      lambda1
      RW
    end

    methods
      function obj = Empleado(nombre, edad, salario, puesto)
        % Llamar al constructor de la superclase
        obj@Persona(nombre, edad)
        obj.Salario = salario;
        obj.Puesto = puesto;
      end
        % Método get para la propiedad dependiente c
        function value1 = get.lambda1(obj)
          value1 = obj.niu;
        end

        function value2 = get.RW(obj)
          value2 = obj.b / obj.NR;
        end

        function mostrarInfo(obj)
          disp(['alfainicio: ' obj.alfainicio])
          disp(['b: ' obj.b])
          disp(['ds: ' obj.ds])
          disp(['dsta: ' obj.dsta])
          disp(['fibdens: ' obj.fibdens])
          disp(['fmc: ' obj.fmc])
          disp(['fvc: ' obj.fvc])
          disp(['niu: ' obj.niu])
          disp(['NP: ' obj.NP])
          disp(['NR: ' obj.NR])
          disp(['resdens: ' obj.resdens])
          disp(['rfm: ' obj.rfm])
          disp(['rim: ' obj.rim])
          disp(['TEX: ' obj.TEX])
          disp(['vel: ' obj.vel])
          disp(['zfg: ' obj.zfg])
          disp(['zfm: ' obj.zfm])
          disp(['zi: ' obj.zi])
          disp(['zig: ' obj.zig])
          disp(['zim: ' obj.zim])
          disp(['lambda1: ' obj.lambda1])
          disp(['RW: ' obj.RW])
        end
        function excelInfo(obj)
        % Escribir una celda, ya que es más fácil exportar un archivo con formato de celdas
        celda_variables_de_inicio = {};
        celda_variables_de_inicio(1,:) = {'alfainicio', obj.alfainicio};
        celda_variables_de_inicio(2,:) = {'b', obj.b};
        celda_variables_de_inicio(3,:) = {'ds', obj.ds};
        celda_variables_de_inicio(4,:) = {'dsta', obj.dsta};
        celda_variables_de_inicio(5,:) = {'fibdens', obj.fibdens};
        celda_variables_de_inicio(6,:) = {'fmc', obj.fmc};
        celda_variables_de_inicio(7,:) = {'fvc', obj.fvc};
        celda_variables_de_inicio(8,:) = {'niu', obj.niu};
        celda_variables_de_inicio(9,:) = {'NP', obj.NP};
        celda_variables_de_inicio(10,:) = {'NR', obj.NR};
        celda_variables_de_inicio(11,:) = {'resdens', obj.resdens};
        celda_variables_de_inicio(12,:) = {'rfm', obj.rfm};
        celda_variables_de_inicio(13,:) = {'rim', obj.rim};
        celda_variables_de_inicio(14,:) = {'TEX', obj.TEX};
        celda_variables_de_inicio(15,:) = {'vel', obj.vel};
        celda_variables_de_inicio(16,:) = {'zfg', obj.zfg};
        celda_variables_de_inicio(17,:) = {'zfm', obj.zfm};
        celda_variables_de_inicio(18,:) = {'zi', obj.zi};
        celda_variables_de_inicio(19,:) = {'zig', obj.zig};
        celda_variables_de_inicio(20,:) = {'zim', obj.zim};
        celda_variables_de_inicio(21,:) = {'lambda1', obj.lambda1};
        celda_variables_de_inicio(22,:) = {'RW', obj.RW};
        
        % Escribir un archivo con las variables de entrada
        writecell(celda_variables_de_inicio, obj.excel, 'Sheet', 'Datos de entrada', 'Range', 'A1');
        end
    end
end