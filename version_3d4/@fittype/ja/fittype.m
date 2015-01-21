%FITTYPE  fit type �I�u�W�F�N�g���쐬
%
%   FITTYPE(EXPR) �́A������ EXPR �̒��Ɋ܂܂�� MATLAB �\������ FITTYPE 
%   �֐��I�u�W�F�N�g���쐬���܂��B���͈����́A�ϐ����ɑ΂��āAEXPR ��
%   ���􂷂邱�ƂŁA�����I�Ɍ��肳��܂� (SYMVAR ���Q��)�B���̏ꍇ�A
%   'x' �͓Ɨ��ϐ��A'y' �͏]���ϐ��ŁA���ׂĂ̑��̕ϐ��̓��f���̌W����
%   ���肵�Ă��܂��B�ϐ������݂��Ȃ��ꍇ�A'x' ���g���܂��B���ׂĂ̌W����
%   �X�J���ł��BEXPR �́A����`���f���Ƃ��Ď�舵���܂��B
%
%   ���̕����� EXPR �Ŏ��̌W�������W���Ƃ��ėp���邱�Ƃ͂ł��܂���B
%   i, j, pi, inf, nan, eps.
%
%   ���`���f���́A�ȉ��̌^�̃��f���ł��B
%      coeff1 * term1 + coeff2 * term2 + coeff3 * term3 + ...
%   (�����ŁA�W���́Aterm1, term2, ���X�̒��ɂ͕\��܂���) ���̃��f���́A
%   �Z���z��Ƃ��āAEXPR ��ݒ肷�邱�ƂŋK�肳��A�W�����܂܂Ȃ����`���f����
%   �e���́AEXPR �̃Z���ɐݒ肳��܂��B���Ƃ��΁A���f�� 'a*x + b*sin(x) + c' �́A
%   'a', 'b', 'c' ���܂񂾌^�Ő��`�ł��B3 �̍��A'x', 'sin(x)', '1' (c=c*1 �̂���) 
%   �������AEXPR �́A3 �̃Z������ɂ܂Ƃ߂��Z�� {'x','sin(x)','1'} �ɂȂ�܂��B
%
%   FITTYPE(TYPE) �́A�^�C�v TYPE �� FITTYPE �I�u�W�F�N�g���쐬���܂��B
%   TYPE �ɑ΂���I���F
%                TYPE                   �ڍ�
%   �X�v���C���F
%                'smoothingspline'      �������X�v���C��
%                'cubicspline'          �O�� (���}���ꂽ) �X�v���C��
%   ���:
%                'linearinterp'         ���`���}
%                'nearestinterp'        �ŋߖT���}
%                'splineinterp'         �O���X�v���C�����}
%                'pchipinterp'          �^��ۑ����� (pchip) ���}
%
%   �܂��́ACFLIBHELP �ɋL�q����郉�C�u�������f���̖��O��ݒ�ł��܂� 
%   (���C�u�������f���̖��O�A�ڍׂ́Atype CFLIBHELP ���Q�Ƃ��Ă�������)�B
%
%   FITTYPE(EXPR,PARAM1,VALUE1,PARAM2,VALUE2,....) �́APARAM-VALUE �̑g��
%   �g���āA�f�t�H���g�l�����������܂��B
%   'independent'    �Ɨ��ϐ������w��
%   'dependent'      �]���ϐ������w��
%   'coefficients'   �W�������w�� (2 �ȏ�̏ꍇ�́A�Z���z��)�B��L��
%                    ���O���ɒ��ӂ��Ă��������B
%   'problem'        ���Ɉˑ�����(�萔)�����w��
%                    (2 �ȏ�̏ꍇ�́A�Z���z��)
%   'options'        ���̕������ɑ΂��ẮA�f�t�H���g�� 'fitoptions' ��ݒ�
%   �f�t�H���g�F �Ɨ��ϐ��� x �ł��B
%                �]���ϐ��́Ay �ł��B
%                ���ɏ]������萔�͂���܂���B
%                ���̑��́A���ׂČW�����ł��B
%
%   �����̕����V���{�������g�����Ƃ��ł��܂��B
%
%   ��:
%     g = fittype('a*x^2+b*x+c')
%     g = fittype('a*x^2+b*x+c','coeff',{'a','b','c'})
%     g = fittype('a*time^2+b*time+c','indep','time')
%     g = fittype('a*time^2+b*time+c','indep','time','depen','height')
%     g = fittype('a*x+n*b','problem','n')
%     g = fittype({'cos(x)','1'})                            % ���`
%     g = fittype({'cos(x)','1'}, 'coefficients', {'a','b'}) % ���`
%
%   �Q�l CFLIBHELP, FIT.


%   Copyright 1999-2008 The MathWorks, Inc.
