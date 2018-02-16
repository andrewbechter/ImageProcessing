function [p] = Centroid(data)

N = size(data,2);

for ii = 1:N
    M = size(data(:,ii),1);

    for jj = 1:M
        frame=data{jj,ii};
        if isempty(frame)==0
            [~,~,x_centroid,y_centroid,sigma_x,sigma_y,max_count,theta] = BasCentroid(frame,M,jj);
            p{jj,ii,1}=x_centroid;
            p{jj,ii,2}=y_centroid;
            p{jj,ii,3}=sigma_x;
            p{jj,ii,4}=sigma_y;
            p{jj,ii,5}=max_count;
            p{jj,ii,6}=theta;
            
        else
            disp('empty')
            p{jj,ii,1}=NaN;
            p{jj,ii,2}=NaN;
            p{jj,ii,3}=NaN;
            p{jj,ii,4}=NaN;
            p{jj,ii,5}=NaN;
            p{jj,ii,6}=NaN;
        end
        
    end
end