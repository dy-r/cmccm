FROM opencpu/base

RUN R -e 'devtools::install_github("dy-r","cmccm")'

RUN \
  echo 'Redirect /index.html /ocpu/library/cmccm/www' > /etc/apache2/sites-enabled/app.conf
