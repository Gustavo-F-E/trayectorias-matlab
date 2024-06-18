classdef variables_de_salida

    properties
        a =9
        b =5
    end

    properties (Dependent)
        c
    end

    methods
        function obj = variables_de_salida(a, b)
            if nargin > 0
                obj.a = a;
                obj.b = b;
            end
        end

        % Método get para la propiedad dependiente c
        function value = get.c(obj)
            value = obj.a + obj.b;
        end

        function mostrarInfo(obj)
          disp(['a: ', num2str(obj.a)])
          disp(['b: ', num2str(obj.b)])
          disp(['c: ', num2str(obj.c)])
      end
    end
end


function vector_variables_de_salida = variables_de_salida(excel)
  vector_variables_de_salida = struct( ...
      'nomeoutput', '04.03.2016_rc75_b10_-500,0,300,1000,alfa20', ...
      'filepath', 'C:\Tesis de maestria\Archivos 0\0-Filament winding\fundamentos de patrones\nuevo_trayectoria_filament_winding', ...
      'gravartrajectorias', 's', ...
      'cnc', 's', ...
      'gravarcnc', 's', ...
      'D', 50, ...
      'seccoes', 'rectangular', ...
      'belip', 1, ...
      'sobre', 'sobreposiçao2');
  % Escribir una celda, ya que es más fácil exportar un archivo con formato de celdas
  celda_variables_de_salida (1,:) = {'nomeoutput',  vector_variables_de_salida.nomeoutput};
  celda_variables_de_salida (2,:) = {'filepath',  vector_variables_de_salida.filepath};         
  celda_variables_de_salida (3,:) = {'gravartrajectorias',  vector_variables_de_salida.gravartrajectorias};        
  celda_variables_de_salida (4,:) = {'cnc',  vector_variables_de_salida.cnc};       
  celda_variables_de_salida (5,:) = {'gravarcnc',  vector_variables_de_salida.gravarcnc};  
  celda_variables_de_salida (6,:) = {'D',  vector_variables_de_salida.D};       
  celda_variables_de_salida (7,:) = {'seccoes',  vector_variables_de_salida.seccoes};        
  celda_variables_de_salida (8,:) = {'belip',  vector_variables_de_salida.belip};      
  celda_variables_de_salida (9,:) = {'sobre',  vector_variables_de_salida.sobre};         
  % Escribir un archivo con las variables de entrada
  writecell(celda_variables_de_salida, excel, 'Sheet', 'Salida de archivos', 'Range', 'A1');
end