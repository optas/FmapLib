classdef Laplacian
    % All the goodies around the Laplacian of a graph.
    % (c) Panos Achlioptas 2014  -  http://www.stanford.edu/~optas/FmapLib
    
    properties (SetAccess = private)
        type      % (String) Describes the type of Laplacian (See constructor for valid values).
        L         % (n x n)  The Laplacian matrix.
        G         % (Graph)  The corresponding graph from which the Laplacian is derived.
    end
            
    methods (Access = public)
        % Class Constructor.
        function obj = Laplacian(varargin)            
            A, laplacian_type
            obj.L      = adjacency_to_laplacian(A, laplacianType);
            obj.type   = laplacianType;                                
        end                   
    end
    
    methods (Static)
        function [L] = adjacency_to_laplacian(A, laplacian_type)    
            % Computes the laplacian matrix for a graph described by its adjacency matrix.
            %  
            % Input:    A                - (n x n) Adjacency matrix of a graph with n nodes.
            %           laplacian_type   - (String, optional) Describes the desired type of laplacian.
            %
            %                           Valid values:                           
            %                               'comb' - Combinatorial (unormalized) laplacian (Default value).
            %                               'norm' - Symmetric Normalized Laplacian.
            %                               'sign' - Signless Laplacian.
            %                                       
            % Output:   L               - (n x n) sparse matrix of the corresponding laplacian.
            %                            
            % Notes:  
            %       DOI: "A Tutorial on Spectral Clustering, U von Luxburg".
            %
            % (c) Panos Achlioptas 2015  -  http://www.stanford.edu/~optas/FmapLib

            if(size(A, 1) ~= size(A, 2))
                error('Adjacency matrices are square, but a non square matrix was given.')
            end

            if nargin < 2
                laplacian_type = 'comb';  % Default output is the combinatorial Laplacian.
            end

            n = size(A, 1);
            total_weight = sum(A, 2);            
            D = spdiags(total_weight, 0, n, n);

            if strcmp(laplacian_type, 'comb')         
                L = -A + D;      
            elseif strcmp(laplacian_type, 'norm')
                Dn = spdiags (1 ./ sqrt(total_weight), 0, n, n);
                L  = Dn *( -A + D) * Dn;
            elseif strcmp(laplacian_type, 'sign')         
                L = A + D;      
            else
                error('Please provide a valid argument for the type of laplacian you want.')
            end
        end
        
    end
       
end

