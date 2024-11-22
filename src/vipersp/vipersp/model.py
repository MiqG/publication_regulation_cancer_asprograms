import torch
import torch.nn as nn
import lightning as L
from torchmetrics.regression import PearsonCorrCoef

class ElementwiseLinear(nn.Module):
    def __init__(self, input_size: int) -> None:
        super(ElementwiseLinear, self).__init__()

        # w is the learnable weight of this layer module
        self.w = nn.Parameter(torch.rand(input_size), requires_grad=True)

    def forward(self, x: torch.tensor) -> torch.tensor:
        # simple elementwise multiplication
        return self.w * x

class EWlayer(nn.Module):
    def __init__(self, input_size):
        super(EWlayer, self).__init__()
        self.ew = ElementwiseLinear(input_size)

    def forward(self, x):
        x = self.ew(x)
        return x

class FClayer(nn.Module):
    def __init__(self, input_size, output_size):
        super(FClayer, self).__init__()
        self.fc = nn.Linear(input_size, output_size)

    def forward(self, x):
        x = self.fc(x)
        return x

## wraper for training with lightning
class LitWrapper(L.LightningModule):
    def __init__(self, model, criterion, learning_rate):
        super(LitWrapper, self).__init__()
        self.model = model
        self.criterion = criterion
        self.learning_rate = learning_rate
        self.losses = []
        self.pearson = PearsonCorrCoef()

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_pred = self(x)
        loss = self.criterion(y_pred, y)
        corr_pearson = self.pearson(y_pred.ravel(), y.ravel())
        self.log('train_loss', loss)
        self.log('train_pearson', corr_pearson)
        self.losses.append(
            {"epoch": self.current_epoch, "batch": batch_idx, "train_loss": float(loss), "train_pearson": float(corr_pearson)}
        )
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_pred = self(x)
        loss = self.criterion(y_pred, y)
        corr_pearson = self.pearson(y_pred.ravel(), y.ravel())
        self.log('val_loss', loss)
        self.log('val_pearson', corr_pearson)
        self.losses.append(
            {"epoch": self.current_epoch, "batch": batch_idx, "val_loss": float(loss), "val_pearson": float(corr_pearson)}
        )
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)