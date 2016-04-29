function [val] = f(states,e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %#1
    temp=diffRow([states(1)+e,states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1)-e,states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val1_1=temp(1);
    val1_2=temp(2);
    val1_3=temp(3);
    val1_4=temp(4);
    val1_5=temp(5);
    val1_6=temp(6);
    val1_7=temp(7);
    val1_8=temp(8);
    val1_9=temp(9);
    val1_10=temp(10);
    
     %#2
    temp=diffRow([states(1),states(2)+e,states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2)-e,states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val2_1=temp(1);
    val2_2=temp(2);
    val2_3=temp(3);
    val2_4=temp(4);
    val2_5=temp(5);
    val2_6=temp(6);
    val2_7=temp(7);
    val2_8=temp(8);
    val2_9=temp(9);
    val2_10=temp(10);
    
        %#3
    temp=diffRow([states(1),states(2),states(3)+e,states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3)-e,states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val3_1=temp(1);
    val3_2=temp(2);
    val3_3=temp(3);
    val3_4=temp(4);
    val3_5=temp(5);
    val3_6=temp(6);
    val3_7=temp(7);
    val3_8=temp(8);
    val3_9=temp(9);
    val3_10=temp(10);
    
        %#4
    temp=diffRow([states(1),states(2),states(3),states(4)+e,states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4)-e,states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val4_1=temp(1);
    val4_2=temp(2);
    val4_3=temp(3);
    val4_4=temp(4);
    val4_5=temp(5);
    val4_6=temp(6);
    val4_7=temp(7);
    val4_8=temp(8);
    val4_9=temp(9);
    val4_10=temp(10);
    
        %#5
    temp=diffRow([states(1),states(2),states(3),states(4),states(5)+e,states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5)-e,states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val5_1=temp(1);
    val5_2=temp(2);
    val5_3=temp(3);
    val5_4=temp(4);
    val5_5=temp(5);
    val5_6=temp(6);
    val5_7=temp(7);
    val5_8=temp(8);
    val5_9=temp(9);
    val5_10=temp(10);
    
        %#6
    temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6)+e,states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6)-e,states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val6_1=temp(1);
    val6_2=temp(2);
    val6_3=temp(3);
    val6_4=temp(4);
    val6_5=temp(5);
    val6_6=temp(6);
    val6_7=temp(7);
    val6_8=temp(8);
    val6_9=temp(9);
    val6_10=temp(10);
    
        %#7
    temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7)+e,states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7)-e,states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val7_1=temp(1);
    val7_2=temp(2);
    val7_3=temp(3);
    val7_4=temp(4);
    val7_5=temp(5);
    val7_6=temp(6);
    val7_7=temp(7);
    val7_8=temp(8);
    val7_9=temp(9);
    val7_10=temp(10);
    
        %#8
    temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8)+e,states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8)-e,states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val8_1=temp(1);
    val8_2=temp(2);
    val8_3=temp(3);
    val8_4=temp(4);
    val8_5=temp(5);
    val8_6=temp(6);
    val8_7=temp(7);
    val8_8=temp(8);
    val8_9=temp(9);
    val8_10=temp(10);
    
        %#9
    temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9)+e,states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9)-e,states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val9_1=temp(1);
    val9_2=temp(2);
    val9_3=temp(3);
    val9_4=temp(4);
    val9_5=temp(5);
    val9_6=temp(6);
    val9_7=temp(7);
    val9_8=temp(8);
    val9_9=temp(9);
    val9_10=temp(10);
    
        %#10
    temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10)+e,...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10)-e,...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val10_1=temp(1);
    val10_2=temp(2);
    val10_3=temp(3);
    val10_4=temp(4);
    val10_5=temp(5);
    val10_6=temp(6);
    val10_7=temp(7);
    val10_8=temp(8);
    val10_9=temp(9);
    val10_10=temp(10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %#11
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11)+e,states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11)-e,states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val11_1=temp(1);
    val11_2=temp(2);
    val11_3=temp(3);
    val11_4=temp(4);
    val11_5=temp(5);
    val11_6=temp(6);
    val11_7=temp(7);
    val11_8=temp(8);
    val11_9=temp(9);
    val11_10=temp(10);
    
    %#12
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12)+e,states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12)-e,states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val12_1=temp(1);
    val12_2=temp(2);
    val12_3=temp(3);
    val12_4=temp(4);
    val12_5=temp(5);
    val12_6=temp(6);
    val12_7=temp(7);
    val12_8=temp(8);
    val12_9=temp(9);
    val12_10=temp(10);
    
    %#13
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13)+e,states(14),states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13)-e,states(14),states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val13_1=temp(1);
    val13_2=temp(2);
    val13_3=temp(3);
    val13_4=temp(4);
    val13_5=temp(5);
    val13_6=temp(6);
    val13_7=temp(7);
    val13_8=temp(8);
    val13_9=temp(9);
    val13_10=temp(10);
    
    %#14
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14)+e,states(15),states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14)-e,states(15),states(16),states(17),states(18),states(19),states(20)],e);
    
    val14_1=temp(1);
    val14_2=temp(2);
    val14_3=temp(3);
    val14_4=temp(4);
    val14_5=temp(5);
    val14_6=temp(6);
    val14_7=temp(7);
    val14_8=temp(8);
    val14_9=temp(9);
    val14_10=temp(10);
    
    %#15
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15)+e,states(16),states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15)-e,states(16),states(17),states(18),states(19),states(20)],e);
    
    val15_1=temp(1);
    val15_2=temp(2);
    val15_3=temp(3);
    val15_4=temp(4);
    val15_5=temp(5);
    val15_6=temp(6);
    val15_7=temp(7);
    val15_8=temp(8);
    val15_9=temp(9);
    val15_10=temp(10);
    
    %#16
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16)+e,states(17),states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16)-e,states(17),states(18),states(19),states(20)],e);
    
    val16_1=temp(1);
    val16_2=temp(2);
    val16_3=temp(3);
    val16_4=temp(4);
    val16_5=temp(5);
    val16_6=temp(6);
    val16_7=temp(7);
    val16_8=temp(8);
    val16_9=temp(9);
    val16_10=temp(10);
    
    %#17
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17)+e,states(18),states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17)-e,states(18),states(19),states(20)],e);
    
    val17_1=temp(1);
    val17_2=temp(2);
    val17_3=temp(3);
    val17_4=temp(4);
    val17_5=temp(5);
    val17_6=temp(6);
    val17_7=temp(7);
    val17_8=temp(8);
    val17_9=temp(9);
    val17_10=temp(10);
    
    %#18
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18)+e,states(19),states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18)-e,states(19),states(20)],e);
    
    val18_1=temp(1);
    val18_2=temp(2);
    val18_3=temp(3);
    val18_4=temp(4);
    val18_5=temp(5);
    val18_6=temp(6);
    val18_7=temp(7);
    val18_8=temp(8);
    val18_9=temp(9);
    val18_10=temp(10);
    
    %#19
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19)+e,states(20)],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19)-e,states(20)],e);
    
    val19_1=temp(1);
    val19_2=temp(2);
    val19_3=temp(3);
    val19_4=temp(4);
    val19_5=temp(5);
    val19_6=temp(6);
    val19_7=temp(7);
    val19_8=temp(8);
    val19_9=temp(9);
    val19_10=temp(10);
    
    %#20
        temp=diffRow([states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                  states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)+e],...
                 [states(1),states(2),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10),...
                 states(11),states(12),states(13),states(14),states(15),states(16),states(17),states(18),states(19),states(20)-e],e);
    
    val20_1=temp(1);
    val20_2=temp(2);
    val20_3=temp(3);
    val20_4=temp(4);
    val20_5=temp(5);
    val20_6=temp(6);
    val20_7=temp(7);
    val20_8=temp(8);
    val20_9=temp(9);
    val20_10=temp(10);
    
    
    
    
    
    val=[[val1_1,val2_1,val3_1,val4_1,val5_1,val6_1,val7_1,val8_1,val9_1,val10_1,val11_1,val12_1,val13_1,val14_1,val15_1,val16_1,val17_1,val18_1,val19_1,val20_1];
         [val1_2,val2_2,val3_2,val4_2,val5_2,val6_2,val7_2,val8_2,val9_2,val10_2,val11_2,val12_2,val13_2,val14_2,val15_2,val16_2,val17_2,val18_2,val19_2,val20_2];
         [val1_3,val2_3,val3_3,val4_3,val5_3,val6_3,val7_3,val8_3,val9_3,val10_3,val11_3,val12_3,val13_3,val14_3,val15_3,val16_3,val17_3,val18_3,val19_3,val20_3];
         [val1_4,val2_4,val3_4,val4_4,val5_4,val6_4,val7_4,val8_4,val9_4,val10_4,val11_4,val12_4,val13_4,val14_4,val15_4,val16_4,val17_4,val18_4,val19_4,val20_4];
         [val1_5,val2_5,val3_5,val4_5,val5_5,val6_5,val7_5,val8_5,val9_5,val10_5,val11_5,val12_5,val13_5,val14_5,val15_5,val16_5,val17_5,val18_5,val19_5,val20_5];
         [val1_6,val2_6,val3_6,val4_6,val5_6,val6_6,val7_6,val8_6,val9_6,val10_6,val11_6,val12_6,val13_6,val14_6,val15_6,val16_6,val17_6,val18_6,val19_6,val20_6];
         [val1_7,val2_7,val3_7,val4_7,val5_7,val6_7,val7_7,val8_7,val9_7,val10_7,val11_7,val12_7,val13_7,val14_7,val15_7,val16_7,val17_7,val18_7,val19_7,val20_7];
         [val1_8,val2_8,val3_8,val4_8,val5_8,val6_8,val7_8,val8_8,val9_8,val10_8,val11_8,val12_8,val13_8,val14_8,val15_8,val16_8,val17_8,val18_8,val19_8,val20_8];
         [val1_9,val2_9,val3_9,val4_9,val5_9,val6_9,val7_9,val8_9,val9_9,val10_9,val11_9,val12_9,val13_9,val14_9,val15_9,val16_9,val17_9,val18_9,val19_9,val20_9];
         [val1_10,val2_10,val3_10,val4_10,val5_10,val6_10,val7_10,val8_10,val9_10,val10_10,val11_10,val12_10,val13_10,val14_10,val15_10,val16_10,val17_10,val18_10,val19_10,val20_10]];
end

