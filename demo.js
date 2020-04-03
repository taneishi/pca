document.addEventListener('DOMContentLoaded', function () {
    var myChart = Highcharts.chart('container', {
        chart: {
            type: 'scatter',
            width: 600,
            height: 600
        },
        title: {
            text: 'Principal Component Anaylsis'
        },
        xAxis: {
            title: {
                text: 'PC1'
            }
        },
        yAxis: {
            title: {
                text: 'PC2'
            }
        },
    });

    d3.csv('iris.csv', function(error, data){
        columns = ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'];
        
        X = data.map(function(row){
            return columns.map(function(col){ return Number(row[col]); })
        });
        
        pc = new PCA().pca(X, 2);

        categories = data.map(function(row){ return row['Species']; });
        pc = pc.map(function(row, i){ return row.concat(categories[i]); });

        categories = d3.set(categories).values();

        categories.map(function(category){
            myChart.addSeries({
                name: category, 
                data: pc.filter(function(row){ return row[4] == category; })
            });
        });
    })
});

