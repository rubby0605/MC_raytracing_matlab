tt=100;
Number=74971;
%% vector of face
aa=zeros(3,1);
bb=aa;cc=aa;c1=aa;c2=aa;c3=aa;
face=zeros(nface,3);
dir=zeros(nface,2);  % location of face (theta,fi)
locxyz=zeros(nface,3);
maxi=zeros(3,1);
maxnum=zeros(nface,3);
face_area=zeros(nface,1);
% Normalize the location of spot <=1
for np=1:npoint, spot(np,:)=spot(np,:)./sqrt(sum(spot(np,:).^2));
end

for nf=1:nface
    
    aa(:)=spot(facet(nf,1),:);
    bb(:)=spot(facet(nf,2),:);
    if(sqrt(sum(aa(:).^2))>=sqrt(sum(bb(:).^2))),maxi(:)=aa(:);
    else maxi(:)=bb(:);
    end
    cc(:)=spot(facet(nf,3),:);
    if(sqrt(sum(maxi(:).^2))<=sqrt(sum(cc(:).^2))),
       maxi(:)=cc(:);
    end
    maxnum(nf,:)=maxi(:);
    locxyz(nf,:)=(aa(:)+bb(:)+cc(:))./3;
    
    c1(:)=bb(:)-aa(:);
    c2(:)=cc(:)-aa(:);
    face_area(nf)=dot(c1,c2)./2;
    c3(:)=cross(c2./sqrt(sum(c2(:).^2)),c1./sqrt(sum(c1(:).^2))); % outward
    face(nf,:)=c3(:)./sqrt(sum(c3(:).^2));
    thi=asin(c3(3)/sqrt(sum(c3(:).^2)));
    if(c3(1)<0),fi=-atan(c3(2)/c3(1))+3*pi/2;
    else fi=-atan(c3(2)/c3(1))+pi/2;
    end
    dir(nf,:)=[thi fi];
    if(sum(c3(:).^2)==0),printf('strange in %d\n',nf);end
end
%face(1:10,:)
%% write the Surface area

str3=['Face_Area_' str '.dat'];
fid=fopen(str3,'w');
fprintf(fid,'%d\n',nface);

for nf=1:nface
 
        fprintf(fid,'%9.7f\n',face_area(nf));
end
fclose(fid);


%% Read the property file
 
 fid=fopen('Property.dat','r');

    line = fgetl(fid);
    total_number=sscanf(line,'%d');
    
    line = fgetl(fid);
    stri=sscanf(line,'%6s %s %s %s %s %s %s');

for Num=1:total_number
    line = fgetl(fid);
    PP=sscanf(line, '%f %f %f %f %f %f %s');
    if(PP(1)==Number),fprintf('The asteroid is %s\n',PP(7:max(size(PP))));break;end
end
fclose(fid);
za=PP(6)
%% Count for the shadow neiborhood
count=zeros(tt,nface,10);
stopcount=zeros(tt,nface);
tic;
initial_sc=0.89;
%copy_sc=zeros(tt,nface);
za=za/180*pi;

for time=1:tt
    fi=(time-0.5)/tt*2*pi;
    solar=[cos(-za)*sin(fi) cos(-za)*cos(fi) sin(-za)];
for nf=1054:1054
    %jj=0;
    %kk=0;
    sc=initial_sc;
    %while(jj<=10)
        jj=0;
        for nnf=1:nface
        io=0;
        if(nnf==nf),continue;end
        s1=(spot(facet(nnf,1),:)+spot(facet(nnf,1),:)+spot(facet(nnf,1),:))/3;
        s2=(spot(facet(nf,1),:)+spot(facet(nf,1),:)+spot(facet(nf,1),:))/3;
        line_s=s1(:)-s2(:);
        line_s=line_s(:)./sqrt(sum(line_s(:).^2));
        if(dot(line_s,solar)>=sc)
            jj=jj+1;
            count(time,nf,jj)=nnf;
        end
        %if(jj>=100),sc=(sc+1)/2;jj=0;break;end
        end
     %   if(jj<100),break;end
     %   disp(kk)
     %   kk=kk+1;
        %if(jj<=),sc=(sc+0)/2;end
   % end
    %disp(sc)
    %copy_sc(time,nf)=sc;
    stopcount(time,nf)=jj;
    disp(nf)
end
disp(time)
end
toc
nf
%% "patch" the Flux & save as figures
TT=zeros(nface,1);
%figure
target_face=60;

for n=1:3
for time=1:tt
    hold off;
    TT=zeros(nface,1);
TT(count(time,target_face,:)>0)=1;
TT(target_face)=0.6;
%TT(:)=CC6(:);
p = patch('Faces',facet,'Vertices',spot,'FaceColor','w');
cdata=TT;
set(p,'FaceColor','flat',...
'CData',cdata,...
'CDataMapping','scaled')
str=['FEM_dust_aa' num2str(time) '.jpg']
colorbar
%saveas(1,str);
%view([-time/2*30 -10]);
axis square
hold on
fi=(time-0.5)/tt*2*pi+pi;
    solar=55.*[cos(za)*sin(fi) cos(za)*cos(fi) sin(za)];
if(n==1),plot3(solar(1),solar(2),solar(3),'rp');end
if(n==2),plot3(solar(1),solar(2),solar(3),'gp');end
if(n==3),plot3(solar(1),solar(2),solar(3),'kp');end
pause(1)

end
hold off;
end

%% Read the Count
depth=10;
str2=['count_dep=' num2str(depth) 'm.dat'];
fid=fopen(str2,'r');
for ii=1:1
    line = fgetl(fid);
    stri=sscanf(line,'%d');
    %str=sscanf(line, '%f');
    if(ii==1),ay=stri;end
end
nface=ay;
count=zeros(nface,10);
stopcount=zeros(nf,1);
for nf=1:nface
    line = fgetl(fid);
    stop=sscanf(line, '%d');
    stopcount(nf)=stop;
    for jj=1:stop
    line = fgetl(fid);
    str=sscanf(line, '%d');
    count(nf,jj)=str;
    end
end
fclose(fid);
%% write it on file
str2=['count_dep=m.dat'];
fid=fopen(str2,'w');
fprintf(fid,'%d\n',nface);
for nf=1:nface
    fprintf(fid,'%d\n',stopcount(nf));
    for fd=1:stopcount(nf)
        fprintf(fid,'%d\n',count(nf,fd));
    end
end
fclose(fid);
whos count


%% shadow
c1=zeros(3,1);
c2=c1;
c3=c1;
c4=c1;
c5=c1;
pot=c1;
da=c1;
db=c1;
dc=c1;
endc=c1;

aa=c1;
bb=c1;
cc=c1;
tt=24;
za=PP(6)
za=za/180*pi;                 % inclination of thita
zb=0; 
F=zeros(nface,tt);
h=0.01;
tic;

%% 
tic;
CC6=zeros(nface,1);
Nmonte=50;
hold on;
for time=1:tt
    fi=(time-0.5)/tt*2*pi+pi;
    solar=[cos(-za)*sin(fi) cos(-za)*cos(fi) sin(-za)];
for nf=1:nface
    c6=face(nf,:);
    c6=dot(solar,c6);
    CC6(nf)=c6;
    if(c6<0),F(nf,time)=0;
    continue;end
    aa(:)=spot(facet(nf,1),:);
    bb(:)=spot(facet(nf,2),:);
    cc(:)=spot(facet(nf,3),:);
    c1(:)=aa(:)-bb(:);
    c2(:)=cc(:)-bb(:);
    cosb=dot(c1./sqrt(sum(c1(:).^2)),c2./sqrt(sum(c2(:).^2)));
    if(cosb<0),copy=aa;aa=bb;bb=copy;c1(:)=aa(:)-bb(:);c2(:)=cc(:)-bb(:);
    cosb=dot(c1./sqrt(sum(c1(:).^2)),c2./sqrt(sum(c2(:).^2)));end
    c3(:)=cosb.*c2(:);
    c4(:)=c1(:)-c3(:);
    c5(:)=c2(:)-c3(:);
    posibility=0;    
    for times=1:Nmonte
        r1=rand;
        r2=rand;
        r3=r1-cosb;
        if(r3>0&&((r2)/(1-r1)>sqrt(sum(c2(:).^2))/(sqrt(sum(c5(:).^2)))))
        pot(:)=(1-r3).*c2(:)+(1-r2).*c4(:)+bb(:);
        elseif(r3<0&&(r2/r1>sqrt(sum(c2(:).^2))/sqrt(sum(c3(:).^2))))
        pot(:)=-r3.*c2(:)+(1-r2).*c4(:)+bb(:); 
        else
        pot(:)=r1.*c2(:)+r2.*c4(:)+bb(:);
        end
        p_io=0;
        if(stopcount(time,nf)==0),continue;end
    for jj=1:stopcount(time,nf)
            nnf=count(time,nf,jj);
            da(:)=spot(facet(nnf,1),:); %point 1 of upper face
            db(:)=spot(facet(nnf,2),:);
            dc(:)=spot(facet(nnf,3),:);
            area=sqrt(sum(cross(db(:)-da(:),dc(:)-da(:)).^2));
            zeropoint=face(nnf,1)*da(1)+face(nnf,2)*da(2)+face(nnf,3)*da(3);
            multi=(zeropoint-(face(nnf,1)*pot(1)+face(nnf,2)*pot(2)+face(nnf,3)*pot(3)))/(face(nnf,1)*solar(1)+face(nnf,2)*solar(2)+face(nnf,3)*solar(3));
            for ijk=1:3
            endc(ijk)=pot(ijk)+solar(ijk).*multi;
            end
            endcda=-endc(:)+da(:);
            endcdb=-endc(:)+db(:);
            endcdc=-endc(:)+dc(:);
            area2=sqrt(sum(cross(endcda,endcdb).^2))+sqrt(sum(cross(endcda,endcdc).^2))+sqrt(sum(cross(endcdc,endcdb).^2));
            if(area2-area<=0.001)
                p_io=1;
                %plot3(endc(1),endc(2),endc(3),'r*');
                %hold on;
                %disp(area2-area);
                %disp(face(nf,1)*endc(1)+face(nf,2)*endc(2)+face(nf,3)*endc(3))
                %if((face(nf,1)*endc(1)+face(nf,2)*endc(2)+face(nf,3)*endc(3))<=0.1),p_io=0;break;end
                break;
            end
     end
        if(p_io==1),p_io=0;posibility=posibility+1;end
     end
    F(nf,time)=(1-posibility/Nmonte).*c6; 
    %disp(F(nf,time))
end
    disp(time)
end

toc

%% write it on file
depth=10;
str3=['Flux_crater.dat'];
fid=fopen(str3,'w');
fprintf(fid,'%d\n',nface);
fprintf(fid,'%d\n',tt);
for nf=1:nface
    for time=1:tt
        fprintf(fid,'%9.7f\n',F(nf,time));
    end
end
fclose(fid);



%% "patch" the Flux & save as figures
TT=zeros(nface,1);
F=abs(F);
for n=1:3
for time=1:tt
    hold off;
%TT(:)=CC6(:);
TT(:)=F(:,time);
p = patch('Faces',facet,'Vertices',spot,'FaceColor','w');
cdata=TT;
set(p,'FaceColor','flat',...
'CData',cdata,...
'CDataMapping','scaled')
str=['FEM_dust_aa' num2str(time) '.jpg']
colorbar
%saveas(1,str);
%view([-time/2*30 -10]);
axis square
hold on
fi=(time-0.5)/tt*2*pi+pi;
    solar=2.*[cos(za)*sin(fi) cos(za)*cos(fi) sin(za)];
if(n==1),plot3(solar(1),solar(2),solar(3),'rp');end
if(n==2),plot3(solar(1),solar(2),solar(3),'gp');end
if(n==3),plot3(solar(1),solar(2),solar(3),'kp');end
pause(0.5)

end
hold off;
end