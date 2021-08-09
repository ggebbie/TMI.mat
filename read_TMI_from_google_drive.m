% Extract data from Google Drive using your favorite method. MATLAB's webread may be an alternative method. Google Drive may ask for spam confirmation for files bigger than 40 MB. Sometimes throws ERROR: cannot verify docs.google.com's certificate, but still works.

% Or, download manually at: https://docs.google.com/uc?export=download&id=1L5i5eQ0QCrrqPKGoAxuB8X4_CltvQBMD and https://docs.google.com/uc?export=download&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm

% Download TMI_4x4deg_data.tar.gz from Google Drive
!wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1L5i5eQ0QCrrqPKGoAxuB8X4_CltvQBMD' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1L5i5eQ0QCrrqPKGoAxuB8X4_CltvQBMD" -O TMI_4x4deg_data.tar.gz && rm -rf /tmp/cookies.txt

% Download TMI_2x2deg_data.tar.gz from Google Drive
! wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm" -O TMI_2x2deg_data.tar.gz && rm -rf /tmp/cookies.txt

!tar xvzf TMI_2x2deg_data.tar.gz 
!rm TMI_2x2deg_data.tar.gz

!tar xvzf TMI_4x4deg_data.tar.gz 
!rm TMI_4x4deg_data.tar.gz
