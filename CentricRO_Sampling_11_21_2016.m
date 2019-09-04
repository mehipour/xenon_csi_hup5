% Name: CentricRO_Sampling_11_21_2016.m
% Purpose: Create a matrix 2 sample a 2D image with an inside-out spiral
% Programmed: 11/17/2016
% Modified:   11/21/2016


% Control Parameters
m_lRepetitionsToMeasure = 15;           % Needs to be incremented by 1 to yield correct loop limit
m_lLinesToMeasure = 16;
lCenterX = uint32((m_lLinesToMeasure+1)/2);
lCenterY = uint32((m_lRepetitionsToMeasure+2)/2);
lRings = uint32((m_lRepetitionsToMeasure+2)/2);
lReorderIndX = 0;
lReorderIndY = 0;
lFilledElemCounter = 0;
lRingElements = 0;
lDir = 0;           % Progress downwards (up = 0, right = 1, down = 2, left = 3)
lIndX = 0;
lIndY = 0;

m_alIndexArrayX  = zeros(m_lRepetitionsToMeasure+1,m_lLinesToMeasure);
m_alIndexArrayY  = zeros(m_lRepetitionsToMeasure+1,m_lLinesToMeasure);
m_alTabXOrdering = zeros(m_lRepetitionsToMeasure+1,m_lLinesToMeasure);
m_alTabYOrdering = zeros(m_lRepetitionsToMeasure+1,m_lLinesToMeasure);

% Re-initialize arrays
m_alIndexArrayX(:, :) = -1;
m_alIndexArrayY(:, :) = -1;

for lI=0:lRings-1
    % Calculate number of elements on each concentric ring
    if (lI == 0)
        lRingElements = 1;
    else
        lRingElements = lI * 8;
    end

    lIndX = lCenterX - lI;
    lIndY = lCenterY;
    
    lDir = 0;        

    for lJ=0:lRingElements-1
        if (lFilledElemCounter < (m_lRepetitionsToMeasure+1) * m_lLinesToMeasure)
            lFilledElemCounter = lFilledElemCounter + 1;
            m_alIndexArrayX(lIndY, lIndX) = lReorderIndX;
            m_alIndexArrayY(lIndY, lIndX) = lReorderIndY;

            if (lDir == 0)       % Up
                if (lIndY == lCenterY - lI)
                    lDir = 1;
                    lIndX = lIndX + 1;
                else
                    lIndY = lIndY - 1;
                end
            elseif (lDir == 1)   % Right
                if (lIndX == lCenterX + lI)
                    lDir = 2;
                    lIndY = lIndY + 1;
                else
                    lIndX = lIndX + 1;
                    if (lIndX > m_lLinesToMeasure)
                        lIndX = 1;
                        lIndY = m_lRepetitionsToMeasure + 1;
                        lDir = 0;
                    end
                end
            elseif (lDir == 2)   % Down
                if (lIndY == lCenterY + lI)
                    lDir = 3;
                    lIndX = lIndX - 1;
                else
                    lIndY = lIndY + 1;
                end
            elseif (lDir == 3)   % Left
                if (lIndX == lCenterX - lI)
                    lDir = 0;
                    lIndY = lIndY - 1;
                else
                    lIndX = lIndX - 1;
                end
            end

            lReorderIndX = lReorderIndX + 1;
            if (lReorderIndX == m_lLinesToMeasure)
                lReorderIndX = lReorderIndX - m_lLinesToMeasure;
                lReorderIndY = lReorderIndY + 1;
            end
        end
    end
end

% Create reordering matrices
for i = 1:m_lRepetitionsToMeasure+1
    for j = 1:m_lLinesToMeasure
        m_alTabXOrdering(m_alIndexArrayY(i,j)+1,m_alIndexArrayX(i,j)+1) = j;
        m_alTabYOrdering(m_alIndexArrayY(i,j)+1,m_alIndexArrayX(i,j)+1) = i;
    end
end
