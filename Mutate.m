function SelCh = Mutate(SelCh,Pm)
%%�������
%����
%SelCh ��ѡ��ĸ���
%Pm �������
%���
%SelCh �����ĸ���

[NSel,L] = size(SelCh);
for i = 1:NSel
    if Pm >rand 
        R = randperm(L);
        if SelCh(i,R(1)) == 0
            SelCh(i,R(1)) = 1;
        else SelCh(i,R(1)) = 0;
        end
    end
end