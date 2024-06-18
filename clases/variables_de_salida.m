classdef variables_de_salida

    properties
      nomeoutput = '04.03.2016_rc75_b10_-500,0,300,1000,alfa20'
      filepath = 'C:\Tesis de maestria\Archivos 0\0-Filament winding\fundamentos de patrones\nuevo_trayectoria_filament_winding'
      gravartrajectorias = 's'
      cnc = 's'
      gravarcnc = 's'
      D = 50
      seccoes = 'rectangular'
      belip = 1
      sobre = 'sobreposiçao2'
      excel ='trayectorias_FW_2.xlsx'
    end

    methods
        function mostrarInfo(obj)
          disp(['nomeoutput: ' obj.nomeoutput])
          disp(['filepath: ' obj.filepath])
          disp(['gravartrajectorias: ' obj.gravartrajectorias])
          disp(['cnc: ' obj.cnc])
          disp(['gravarcnc: ' obj.gravarcnc])
          disp(['D: ' num2str(obj.D)])
          disp(['seccoes: ' obj.seccoes])
          disp(['belip: ' num2str(obj.belip)])
          disp(['sobre: ' obj.sobre])
          disp(['excel: ' obj.excel])
        end
        function excelInfo(obj)
        % Escribir una celda, ya que es más fácil exportar un archivo con formato de celdas
        celda_variables_de_salida = {};
        celda_variables_de_salida (1,:) = {'nomeoutput',  obj.nomeoutput};
        celda_variables_de_salida (2,:) = {'filepath',  obj.filepath};         
        celda_variables_de_salida (3,:) = {'gravartrajectorias',  obj.gravartrajectorias};        
        celda_variables_de_salida (4,:) = {'cnc',  obj.cnc};       
        celda_variables_de_salida (5,:) = {'gravarcnc',  obj.gravarcnc};  
        celda_variables_de_salida (6,:) = {'D',  obj.D};       
        celda_variables_de_salida (7,:) = {'seccoes',  obj.seccoes};        
        celda_variables_de_salida (8,:) = {'belip',  obj.belip};      
        celda_variables_de_salida (9,:) = {'sobre',  obj.sobre};         
        % Escribir un archivo con las variables de entrada
        writecell(celda_variables_de_salida, obj.excel, 'Sheet', 'Salida de archivos', 'Range', 'A1');
        end
    end
end
