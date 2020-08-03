
[img1P]=textread('img1P.txt'); % extracting selected points
[img2P]=textread('img2P.txt');
 
img1=imread('img1.png'); % reading images
img2=imread('img2.png');
 
figure; % ploting images with selected points
subplot(1,2,1),imshow(img1);hold on;
subplot(1,2,1),plot(img1P(:,1),img1P(:,2),'r+');
subplot(1,2,2),imshow(img2);hold on;
subplot(1,2,2),plot(img2P(:,1),img2P(:,2),'r+');
 
%% Calculating F matrix using 8 point algorithm
u1=img1P(:,1);v1=img1P(:,2);
u2=img2P(:,1);v2=img2P(:,2);
 
o=ones(size(u1));
 
A=[u1.*u2 u1.*v2 u1 v1.*u2 v1.*v2 v1 u2 v2 o];
[~,~,V]=svd(A,0);
v=V(:,end);
Pt8F=reshape(v,3,3);
 
 
%% Calculating Normalized values of Points
centroid1 = sum(img1P(:,:))/24;
diff1 = img1P(:,:) - repmat(centroid1,24,1);
avgdist1 = sum(sqrt(diff1(:,1).^2 + diff1( :,2).^2))/24;
scale1 = sqrt(2)/avgdist1;
T1 = diag([scale1 scale1 1]) * [eye(3,2) [-centroid1 1]'];
normpts1 = img1P;
normpts1(:,3) = 1;
normpts1 = normpts1 * T1';
normpts1 = normpts1 ./ repmat(normpts1(:,3), 1, 3);
img1PN=normpts1(:,1:2);
 
centroid2 = sum(img2P(:,:)) /24;
diff2 = img2P(:,:) - repmat(centroid2,24,1);
avgdist2 = sum(sqrt(diff2(:,1).^2 + diff2( :,2).^2))/24;
scale2 = sqrt(2)/avgdist2;
T2 = diag([scale2 scale2 1]) * [eye(3,2) [-centroid2 1]'];
normpts2 = img2P;
normpts2(:,3) = 1;
normpts2 = normpts2 * T2';
normpts2 = normpts2 ./ repmat(normpts2(:,3), 1, 3);
img2PN=normpts2(:,1:2);
 
%% Calculating F matrix using Normalized 8 Point Algorithm
u1=img1PN(:,1); v1=img1PN(:,2);
u2=img2PN(:,1); v2=img2PN(:,2);
o=ones(size(u1));
 
A2=[u1.*u2 u1.*v2 u1 v1.*u2 v1.*v2 v1 u2 v2 o];
[~,~,V]=svd(A2,0);
F = reshape( V(:,end),3 ,3);
[U,S,V] = svd(F);
S(3,3) = 0;             
S = S / S(1,1);         
F = U * S * V'; 
NormPt8F = T2' * F * T1;
 
%% Calculating and Displaying epipolar lines
figure;
imshow(img2);  
hold on;
plot(img2P(:,1),img2P(:,2),'r+');
hold on;
for i=1:1:24
    x=0:1:640;
    L = NormPt8F*[img1P(i,:) 1]';
    y = (-L(3) - L(1)*x) / L(2);
    plot(x,y,'b-');
    hold on;
end
figure;
imshow(img1);
hold on;
plot(img1P(:,1),img1P(:,2),'r+');
hold on;
for i=1:1:24
    x=0:1:640;
    L = NormPt8F*[img2P(i,:) 1]';
    y = (-L(3) - L(1)*x) / L(2);
    plot(x,y,'b-');
    hold on;
end
 
Error_Pt8F = norm( A * reshape(Pt8F',9,1 ));
Error_NormPt8F = norm( A2 * reshape(NormPt8F',9,1 ));
 
%% Calling to Display on Main Window
Pt8F
NormPt8F
Error_Pt8F
Error_NormPt8F
%% End
