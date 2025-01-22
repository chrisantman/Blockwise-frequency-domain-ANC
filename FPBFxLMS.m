%% Algoritmo FxLMS Multicanal Frecuencia, por bloques y particionado
% Sólo permite una señal de ruido por nodo
function [node]=FPBFxLMS(node)

    [IN,B,estado,IN_e]=get_dataIn(node);

    % Leer campos constantes del algortimo
    P=estado.P;
    F=estado.F;
    tam_fft=2*B;
    mu=estado.mu;
    W=estado.W;
    C=estado.C;
    mic = get_micNode(node);
    alt = get_altNode(node);
    K_nodo=length(mic);
    J_nodo=length(alt);
    % Leer campos variables del algortimo
    buffx=estado.buffx;
    x_p=estado.x_p;
    V=estado.V;
    phi=estado.phi;
    x_f=estado.x_f;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        % Buffer de tamaño bloque de la señal de entrada   
        buffx=[buffx(B+1:tam_fft);IN];

        % Filtrado adaptativo
        x_p(:,2:end)=x_p(:,1:P-1);
        x_p(:,1)=buffx;  
        Xp=fft(x_p,tam_fft);

        Y=zeros(tam_fft,P,J_nodo);
        for j=1:J_nodo      	
            Y(:,:,j)=Y(:,:,j)+Xp.*W(:,:,j);
        end
        y=ifft(Y); 
        y=sum(y,2);
        OUT=y(B+1:tam_fft,:);

        % Señal captada por los micros
        IN_e=permute(IN_e,[1 3 2]);
        E=fft([zeros(B,1,K_nodo);IN_e],tam_fft);   
        Ep=repmat(E,1,P); 

       % ADAPTAMOS COEFICIENTES
        for j=1:J_nodo
          v_aux=zeros(tam_fft,P);
          v=V(:,:,1+K_nodo*(j-1):K_nodo*j);
          for k=1:K_nodo 
              v_aux=v_aux+Ep(:,:,k).*conj(v(:,:,k));
          end
          step=ifft(v_aux);
          O=fft([step(1:B,:);zeros(B,P)],tam_fft);
          W(:,:,j)=W(:,:,j)-2*mu*O;
        end

        % Filtrado estima camino secundario
        x_f(:,2:end)=x_f(:,1:F-1);
        x_f(:,1)=buffx;
        Xf=fft(x_f,tam_fft);
        V_aux=zeros(tam_fft,F,K_nodo*J_nodo); 
        for k=1:K_nodo
           for j=1:J_nodo
             V_aux(:,:,k+K_nodo*(j-1))=Xf.*C(:,:,k+K_nodo*(j-1));    
           end
        end  
        Vaux=sum(V_aux,2);   
        V(:,2:end,:)=V(:,1:P-1,:);
        V(:,1,:)=Vaux;        


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Almacenar estados del nodo en estructura
    estado.buffx=buffx;
    estado.x_p=x_p;
    estado.x_f=x_f;
    estado.phi=phi;
    estado.V=V;
    estado.W=W;

    node=set_dataOut(node,OUT,estado);

